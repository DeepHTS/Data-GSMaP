import os
from ftplib import FTP
import datetime
from tqdm import tqdm
import urllib

import numpy as np
from osgeo import gdal, osr
import rasterio
from rasterio.mask import mask
from shapely.geometry import box
from fiona.crs import from_epsg
import geopandas as gpd
import pycrs

from src.helper import argwrapper, imap_unordered_bar, transfer_to_s3

DIR_PARENT_RAW_LOCAL = 'data/GSMaP/raw'
DIR_PARENT_CONVERTED_LOCAL = 'data/GSMaP/converted'
DIR_PARENT_PICKED = 'data/GSMaP/picked'

FTP_ADDRESS_JAXA = 'ftp.gportal.jaxa.jp'
HDF_EXT = ['.h5', '.hdf', '.HDF']
ESPG_GRID = 4326
DICT_GDT = {
    'uint16': gdal.GDT_UInt16,
    'uint8': gdal.GDT_Byte,
    'complex64': gdal.GDT_CFloat64,
    'float32': gdal.GDT_Float32,
    'float64': gdal.GDT_Float64,
    'int16': gdal.GDT_Int16,
    'int32': gdal.GDT_Int32,
    'uint32': gdal.GDT_UInt32,
}

JAXA_FTP_LAYER = ['category', 'project', 'sensor', 'product', 'version']

CATEGORY_STANDARD_GSMAP = 'standard'
CATEGORY_REALTIME_GSMAP = 'nrt'
PROJECT_GSMAP = 'GSMaP'
SENSOR_HOURLY_GSMAP = '3.GSMAP.H'
SENSOR_MONTHLY_GSMAP = '3.GSMAP.M'


def init(*credentials):
    """ set variable of FTP connection as global variable for multiprocessing

    Args:
        *credentials (): [ftp address, user name, password] (for JAXA, 'anonymous' is OK)

    Returns:

    """
    global ftp
    server, user, password = credentials
    ftp = FTP(server)
    ftp.login(user=user, passwd=password)


def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]


def shapes_from_bbox(x_min, y_min, x_max, y_max, epsg_code=4326):
    bbox = box(x_min, y_min, x_max, y_max)
    geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=from_epsg(epsg_code))
    coords = getFeatures(geo)
    return coords


def pick_part_raster(path_raster, path_raster_out, shapes):
    src = rasterio.open(path_raster)
    out_img, out_transform = mask(src, shapes=shapes, crop=True)
    out_meta = src.meta.copy()
    epsg_code = int(src.crs.data['init'][5:])
    out_meta.update({"driver": "GTiff",
                     "height": out_img.shape[1],
                     "width": out_img.shape[2],
                     "transform": out_transform,
                     "crs": pycrs.parse.from_epsg_code(epsg_code).to_proj4()})
    with rasterio.open(path_raster_out, "w", **out_meta) as dest:
        dest.write(out_img)

    return


class DataManagerJAXABase(object):
    def __init__(self,
                 user,
                 dir_parent_raw_local=DIR_PARENT_RAW_LOCAL,
                 dir_parent_converted_local=DIR_PARENT_CONVERTED_LOCAL,
                 processes=1):
        self.user = user
        self.processes = processes
        self.dir_parent_raw_local = dir_parent_raw_local
        self.dir_parent_converted_local = dir_parent_converted_local

        self.ftp_address = FTP_ADDRESS_JAXA
        self.ext = HDF_EXT
        self.dict_gdt = DICT_GDT

    # def __init__(self,
    #              dir_local_parent=DIR_LOCAL_PARENT_ORIGINAL,
    #              save_s3=False,
    #              dir_s3_parent=DIR_S3_PARENT_ORIGINAL,
    #              remove_local_files=True,
    #              dict_ext_filter=DICT_EXT_FILTER,
    #              processes=1):

    def download_from_ftp(self, path_ftp='/standard/GSMaP/3.GSMAP.M/03F/2016/GPMMRG_MAP_1601_M_L3S_MCM_03F.h5'):
        path_local = os.path.join(self.dir_parent_raw_local, path_ftp)

        # make dir
        dir_local = os.path.dirname(path_local)
        if not os.path.exists(dir_local):
            os.makedirs(dir_local)

        # # download from FTP to local
        # if self.processes == 1:
        #     ftp = FTP(self.ftp_address, user=self.user, passwd='anonymous')

        try:
            with open(path_local, 'wb') as f:
                ftp.retrbinary('RETR {0}'.format(path_ftp), f.write)
        except (urllib.error.URLError, FileNotFoundError) as e:
            print('FAILED: {0}'.format(path_ftp))

        return path_local

    def _check_layer(self, band, list_band_filter):
        layer_band = [a for a in band.split('/') if a != '']

        for filter_single in list_band_filter:
            layer_filters = filter_single.split('/')
            layer_filters = [a for a in layer_filters if a != '']
            count = [1 if layer_band[i] != layer_filter else 0 for i, layer_filter in enumerate(layer_filters)]
            if sum(count) == 0:
                return True
        return False

    def _get_filtered_list(self, list_band, list_band_filter):
        if list_band_filter is None:
            return [*range(len(list_band))]
        out = []
        for index, target in enumerate(list_band):
            bool_out = self._check_layer(target, list_band_filter)
            if bool_out:
                out.append(index)
        return out

    def _conv_band_grid(self, band_ds, dir_converted_local, filename_head, nodata=None):
        espg_grid = ESPG_GRID

        band_array = band_ds.ReadAsArray()
        data_type = band_array.dtype.name

        band_array = np.flipud(band_array.T)

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(espg_grid)

        # projection data
        tup_transform = (-180.0, 0.1, 0.0, 90.0, 0.0, -0.1)

        # set filename and path
        filename_grid = filename_head + '_' + str(espg_grid) + '.tif'
        path_grid = os.path.join(dir_converted_local, filename_grid)

        # convert to sinusoidal and save it
        out_ds = gdal.GetDriverByName('GTiff').Create(path_grid,
                                                      band_ds.RasterYSize,
                                                      band_ds.RasterXSize,
                                                      1,  # Number of bands
                                                      self.dict_gdt.get(data_type, gdal.GDT_Unknown),
                                                      ['TILED=YES'])
        out_ds.SetGeoTransform(tup_transform)
        out_ds.SetProjection(srs.ExportToWkt())
        out_ds.GetRasterBand(1).WriteArray(band_array)
        out_ds.FlushCache()
        out_ds = None

        with rasterio.open(path_grid, 'r+') as src:
            src.nodata = nodata

        return path_grid

    def _get_meta_data(self, dataset, dict_meta):
        pass

    def convert_from_hdf_to_gtiff(self, path_local, list_band_filter=['//Grid/monthlyPrecipRate']):
        # strike out not hdf file
        if os.path.splitext(path_local)[-1] not in self.ext:
            print("no HDF or H5 file \n {}".format(path_local))
            dict_path_list = {
                'filepath': path_local,
            }
            return dict_path_list

        dir_converted_local = os.path.join(self.dir_parent_converted_local,
                                           os.path.dirname(path_local).split(self.dir_parent_raw_local)[-1][1:])

        hdf_ds = gdal.Open(path_local)
        try:
            sub_datasets = hdf_ds.GetSubDatasets()
        except AttributeError:
            print("can't open sub datasets by GDAL \n {}".format(path_local))
            return

        # get all layers
        list_band = [dataset[0].split(':')[-1] for dataset in sub_datasets]
        list_index = self._get_filtered_list(list_band, list_band_filter)
        # print([list_band[index] for index in list_index])

        # check exist dest dir
        if not os.path.exists(dir_converted_local):
            os.makedirs(dir_converted_local)

        filename_head = os.path.splitext(os.path.basename(path_local))[0]
        # convert HDF to local
        dict_meta = hdf_ds.GetMetadata()

        list_path_converted = []
        for index in list_index:
            dataset = sub_datasets[index]
            band_ds = gdal.Open(dataset[0], gdal.GA_ReadOnly)
            dict_meta_sub = self._get_meta_data(dataset, dict_meta)
            head = dict_meta_sub.get('head', None)

            if head is not None:
                filename_head_temp = filename_head + '_' + head
            else:
                filename_head_temp = str(filename_head)

            if dict_meta_sub['type'] in ['Grid']:
                path_converted = self._conv_band_grid(band_ds, dir_converted_local, filename_head=filename_head_temp,
                                                      nodata=dict_meta_sub['nodata'])
                list_path_converted.append(path_converted)
            else:
                print("not type: ", dict_meta_sub['type'])

        return list_path_converted


class DataManagerJAXAGSMaP(DataManagerJAXABase):
    def __init__(self, user, dir_parent_raw_local=DIR_PARENT_RAW_LOCAL,
                 dir_parent_converted_local=DIR_PARENT_CONVERTED_LOCAL, processes=1):
        super().__init__(user, dir_parent_raw_local, dir_parent_converted_local, processes)

        self.category_standard = CATEGORY_STANDARD_GSMAP
        self.category_realtime = CATEGORY_REALTIME_GSMAP
        self.project = PROJECT_GSMAP
        self.sensor_hourly = SENSOR_HOURLY_GSMAP
        self.sensor_monthly = SENSOR_MONTHLY_GSMAP

    def _get_meta_data(self, dataset, dict_meta):
        head = '_'.join(dataset[0].split('://')[1].split('/'))
        key_candidate = ['_FillValue', 'CodeMissingValue']
        key = [key for key in dict_meta.keys() if key in [head + '_' + key for key in key_candidate]][0]
        nodata = float(dict_meta[key])
        dict_out = {'nodata': nodata,
                    'head': head,
                    'type': dataset[0].split('://')[1].split('/')[0]}
        return dict_out

    def _filter_filename(self, filename, start_datetime=None, end_datetime=None, ext=None):
        if (ext is not None) and (os.path.splitext(filename)[-1] not in ext):
            return False

        if start_datetime is None:
            start_datetime = datetime.datetime(year=1983, month=1, day=5)
        if end_datetime is None:
            end_datetime = datetime.datetime(year=2100, month=1, day=5)

        product_flag = filename.split('_')[3]
        datetime_flag = filename.split('_')[2]
        if product_flag == 'H':
            year = int('20' + datetime_flag[:2])
            month = int(datetime_flag[2:4])
            day = int(datetime_flag[4:6])
            hour = int(datetime_flag[6:8])
            file_datetime = datetime.datetime(year=year, month=month, day=day, hour=hour)
        elif product_flag == 'M':
            start_datetime = datetime.datetime(year=start_datetime.year, month=start_datetime.month, day=1)
            end_datetime = datetime.datetime(year=end_datetime.year, month=end_datetime.month, day=1)
            year = int('20' + datetime_flag[:2])
            month = int(datetime_flag[2:4])
            file_datetime = datetime.datetime(year=year, month=month, day=1)
        else:
            print('error')
            return

        return file_datetime >= start_datetime and file_datetime <= end_datetime

    def _get_ftp_base_dir(self, ftp_base_dir, category, product, start_datetime=None, end_datetime=None):
        # if 'ftp' not in globals():
        #     ftp = FTP(self.ftp_address, user=self.user, passwd='anonymous')

        if start_datetime is None:
            start_datetime = datetime.datetime(year=1983, month=1, day=5)
        if end_datetime is None:
            end_datetime = datetime.datetime(year=2100, month=1, day=5)

        list_ftp_base_dir = []
        if category == 'standard':
            if product == 'hourly':
                list_ftp_base_dir = ftp.nlst(ftp_base_dir)
                list_ftp_base_dir = [ftp_base_dir for ftp_base_dir in list_ftp_base_dir
                                     if int(os.path.basename(ftp_base_dir)) >= start_datetime.year and
                                     int(os.path.basename(ftp_base_dir)) >= start_datetime.year <= end_datetime.year]
                # month
                list_ftp_base_dir_out = []
                start_year_month = datetime.datetime(year=start_datetime.year, month=start_datetime.month, day=1)
                end_year_month = datetime.datetime(year=end_datetime.year, month=end_datetime.month, day=1)
                for ftp_base_dir in list_ftp_base_dir:
                    year = int(os.path.basename(ftp_base_dir))
                    list_dir_temp = ftp.nlst(ftp_base_dir)
                    for dir_temp in list_dir_temp:
                        month = int(os.path.basename(dir_temp))
                        date_temp = datetime.datetime(year=year, month=month, day=1)
                        if date_temp >= start_year_month and date_temp <= end_year_month:
                            list_ftp_base_dir_out.append(dir_temp)
                list_ftp_base_dir = list(list_ftp_base_dir_out)
                # day
                list_ftp_base_dir_out = []
                for ftp_base_dir in list_ftp_base_dir:
                    month = int(os.path.basename(ftp_base_dir))
                    year = int(ftp_base_dir.split('/')[-2])
                    list_dir_temp = ftp.nlst(ftp_base_dir)
                    for dir_temp in list_dir_temp:
                        day = int(os.path.basename(dir_temp))
                        date_temp = datetime.datetime(year=year, month=month, day=day)
                        if date_temp >= start_datetime and date_temp <= end_datetime:
                            list_ftp_base_dir_out.append(dir_temp)
                list_ftp_base_dir = list(list_ftp_base_dir_out)

            elif product == 'monthly':
                list_ftp_base_dir = ftp.nlst(ftp_base_dir)
                list_ftp_base_dir = [ftp_base_dir for ftp_base_dir in list_ftp_base_dir
                                     if int(os.path.basename(ftp_base_dir)) >= start_datetime.year and
                                     int(os.path.basename(ftp_base_dir)) >= start_datetime.year <= end_datetime.year]
        elif category == 'realtime':
            list_ftp_base_dir = ftp.nlst(ftp_base_dir)

        return list_ftp_base_dir

    def get_path_filtered(self, ftp_dir, start_datetime=None, end_datetime=None, ext=None):
        list_ftp_path_temp = ftp.nlst(ftp_dir)
        list_ftp_path_out = [ftp_path for ftp_path in list_ftp_path_temp
                             if self._filter_filename(os.path.basename(ftp_path),
                                                      start_datetime=start_datetime,
                                                      end_datetime=end_datetime,
                                                      ext=ext)]
        return list_ftp_path_out

    def get_ftp_path_list(self, category='standard', product='hourly', list_version=None, start_datetime=None,
                          end_datetime=None, ext=HDF_EXT):
        # download from FTP to local
        # if 'ftp' not in globals():
        #     ftp = FTP(self.ftp_address, user=self.user, passwd='anonymous')
        credentials = [self.ftp_address, self.user, 'anonymous']
        init(*credentials)

        if category == 'standard':
            path_base = os.path.join(self.category_standard, self.project)
            if product == 'hourly':
                path_base = os.path.join(path_base, self.sensor_hourly)
            elif product == 'monthly':
                path_base = os.path.join(path_base, self.sensor_monthly)
            else:
                print('check args product')
                return

            # get versions:
            if list_version is None:
                list_ftp_base_dir = ftp.nlst(path_base)
            else:
                list_ftp_base_dir = [os.path.join(path_base, version) for version in list_version]

        elif category == 'realtime':
            path_base = os.path.join(self.category_realtime, self.project, self.sensor_hourly)
            list_ftp_base_dir = [path_base]
        else:
            print('check args category')
            return

        # get ftp list under designated criteria
        if self.processes == 1:
            list_ftp_dir = []
            for ftp_base_dir in tqdm(list_ftp_base_dir, total=len(list_ftp_base_dir)):
                list_ftp_dir.extend(
                    self._get_ftp_base_dir(ftp_base_dir, category, product, start_datetime, end_datetime))
        else:
            # credentials = [self.ftp_address, self.user, 'anonymous']
            # init(*credentials)
            func_args = [(self._get_ftp_base_dir, ftp_base_dir, category, product, start_datetime, end_datetime)
                         for ftp_base_dir in list_ftp_base_dir]
            list_ftp_dir = imap_unordered_bar(argwrapper, func_args, self.processes, extend=True, init=init,
                                              credentials=credentials)

        if len(list_ftp_base_dir) == 0:
            return

        # get ftp path list
        if self.processes == 1:
            list_ftp_path = []
            for ftp_dir in tqdm(list_ftp_dir, total=len(list_ftp_dir)):
                list_ftp_path.extend(self.get_path_filtered(ftp_dir, start_datetime=start_datetime,
                                                            end_datetime=end_datetime,
                                                            ext=ext))
        else:
            # credentials = [self.ftp_address, self.user, 'anonymous']
            # init(*credentials)
            func_args = [(self.get_path_filtered, ftp_dir, start_datetime, end_datetime, ext)
                         for ftp_dir in list_ftp_dir]
            list_ftp_path = imap_unordered_bar(argwrapper, func_args, self.processes, extend=True, init=init,
                                               credentials=credentials)

        return list_ftp_path

    def _exec_get_raster(self, ftp_path, convert_to_gtiff=True, list_band_filter=['//Grid/monthlyPrecipRate'],
                         pick_part=True,
                         x_min=120.61, y_min=22.29, x_max=151.35, y_max=46.8, epsg_code=4326, area_name='japan',
                         dir_parent_picked_local=DIR_PARENT_PICKED,
                         save_s3=False, dir_s3_parent=None, remove_local_files=False, s3_bucket_name=None):

        list_out_path = []
        path_local = self.download_from_ftp(ftp_path)

        if not os.path.exists(dir_parent_picked_local):
            os.makedirs(dir_parent_picked_local)

        if convert_to_gtiff:
            list_path_converted = self.convert_from_hdf_to_gtiff(path_local, list_band_filter=list_band_filter)
            if pick_part:
                shapes = shapes_from_bbox(x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max, epsg_code=epsg_code)
                for path_converted in list_path_converted:
                    filename = os.path.basename(path_converted)
                    filename = os.path.splitext(filename)[0] + '_' + area_name + os.path.splitext(filename)[1]
                    path_converted_out = os.path.join(dir_parent_picked_local, filename)
                    pick_part_raster(path_converted, path_converted_out, shapes)
                    list_out_path.append(path_converted_out)
        else:
            list_path_converted = []

        if save_s3:
            if dir_s3_parent is None:
                dir_parent_raw_s3 = str(self.dir_parent_raw_local)
                dir_parent_converted_s3 = str(self.dir_parent_converted_local)
                dir_parent_picked_s3 = str(dir_parent_picked_local)
            else:
                dir_parent_raw_s3 = os.path.join(dir_s3_parent, self.dir_parent_raw_local)
                dir_parent_converted_s3 = os.path.join(dir_s3_parent, self.dir_parent_converted_local)
                dir_parent_picked_s3 = os.path.join(dir_s3_parent, dir_parent_picked_local)

            transfer_to_s3(path_local, dir_local_parent=self.dir_parent_raw_local,
                           dir_s3_parent=dir_parent_raw_s3,
                           remove_local_file=remove_local_files,
                           multiprocessing=self.processes > 1, s3_bucket_name=s3_bucket_name)
            if len(list_path_converted) > 0:
                for path_converted in list_path_converted:
                    transfer_to_s3(path_converted, dir_local_parent=self.dir_parent_converted_local,
                                   dir_s3_parent=dir_parent_converted_s3,
                                   remove_local_file=remove_local_files,
                                   multiprocessing=self.processes > 1, s3_bucket_name=s3_bucket_name)
            if len(list_out_path) > 0:
                for path in list_out_path:
                    transfer_to_s3(path, dir_local_parent=dir_parent_picked_local, dir_s3_parent=dir_parent_picked_s3,
                                   remove_local_file=remove_local_files,
                                   multiprocessing=self.processes > 1, s3_bucket_name=s3_bucket_name)

        return []

    def _makedirs_from_ftp_path(self, list_ftp_path, convert_to_gtiff=True):
        list_dirs = []
        for path_ftp in list_ftp_path:
            path_local = os.path.join(self.dir_parent_raw_local, path_ftp)
            list_dirs.append(os.path.dirname(path_local))

            if convert_to_gtiff:
                dir_converted_local = os.path.join(self.dir_parent_converted_local,
                                                   os.path.dirname(path_local).split(self.dir_parent_raw_local)[-1][1:])
                list_dirs.append(dir_converted_local)
        list_dirs_unique = list(np.unique(list_dirs))

        for dir in list_dirs_unique:
            if not os.path.exists(dir):
                os.makedirs(dir)

        return

    def get_raster_data(self, category='standard', product='hourly', list_version=None, start_datetime=None,
                        end_datetime=None, ext=HDF_EXT, convert_to_gtiff=True,
                        list_band_filter=['//Grid/monthlyPrecipRate'],
                        pick_part=True,
                        x_min=120.61, y_min=22.29, x_max=151.35, y_max=46.8, epsg_code=4326, area_name='japan',
                        dir_parent_picked_local=DIR_PARENT_PICKED,
                        save_s3=False, dir_s3_parent=None, remove_local_files=False, s3_bucket_name=None
                        ):
        list_ftp_path = self.get_ftp_path_list(category=category, product=product, list_version=list_version,
                                               start_datetime=start_datetime, end_datetime=end_datetime, ext=ext)

        self._makedirs_from_ftp_path(list_ftp_path, convert_to_gtiff=convert_to_gtiff)

        # don't need to touch other coords
        # credentials = [self.ftp_address, self.user, 'anonymous']
        # init(*credentials)

        dir_parent_picked_local = os.path.join(dir_parent_picked_local, category, product)
        if pick_part and not os.path.exists(dir_parent_picked_local):
            os.makedirs(dir_parent_picked_local)

        if self.processes == 1:
            list_out_path = []
            # global ftp
            # ftp = FTP(self.ftp_address, user=self.user, passwd='anonymous')
            for ftp_path in tqdm(list_ftp_path, total=len(list_ftp_path)):
                self._exec_get_raster(ftp_path, convert_to_gtiff, list_band_filter, pick_part,
                                      x_min, y_min, x_max, y_max, epsg_code, area_name,
                                      dir_parent_picked_local, save_s3, dir_s3_parent,
                                      remove_local_files, s3_bucket_name)
        else:
            credentials = [self.ftp_address, self.user, 'anonymous']
            init(*credentials)
            func_args = [(self._exec_get_raster, ftp_path, convert_to_gtiff, list_band_filter, pick_part,
                          x_min, y_min, x_max, y_max, epsg_code, area_name,
                          dir_parent_picked_local, save_s3, dir_s3_parent,
                          remove_local_files, s3_bucket_name)
                         for ftp_path in list_ftp_path]
            imap_unordered_bar(argwrapper, func_args, self.processes, extend=True, init=init, credentials=credentials)
        return
