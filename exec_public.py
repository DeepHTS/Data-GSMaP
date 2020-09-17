import datetime
import os

from src.download import DataManagerJAXAGSMaP, shapes_from_bbox

JAXA_GPORTAL_USER = 'xxxx'
S3_BUCKET_NAME = 'xxxx'

def main():
    processes = 12
    start_datetime = datetime.datetime(year=2020, month=4, day=25)
    end_datetime = datetime.datetime(year=2020, month=4, day=26)
    category = 'standard'
    product = 'hourly'
    convert_to_gtiff = True
    list_band_filter = ['//Grid/hourlyPrecipRate']
    pick_part = True
    x_min = 120.61
    y_min = 22.29
    x_max = 151.35
    y_max = 46.8
    epsg_code = 4326
    area_name = 'japan'
    save_s3 = False
    remove_local_files = False

    data_manager_jaxa_gsmap = DataManagerJAXAGSMaP(user=JAXA_GPORTAL_USER, processes=processes)
    list_out_path = data_manager_jaxa_gsmap.get_raster_data(category=category,
                                                            product=product,
                                                            start_datetime=start_datetime,
                                                            end_datetime=end_datetime,
                                                            convert_to_gtiff=convert_to_gtiff,
                                                            list_band_filter=list_band_filter,
                                                            pick_part=pick_part,
                                                            x_min=x_min,
                                                            y_min=y_min,
                                                            x_max=x_max,
                                                            y_max=y_max,
                                                            epsg_code=epsg_code,
                                                            area_name=area_name,
                                                            save_s3=save_s3,
                                                            remove_local_files=remove_local_files,
                                                            s3_bucket_name=S3_BUCKET_NAME)

    return

if __name__ == "__main__":
    main()

    # data_manager_jaxa_gsmap = DataManagerJAXAGSMaP(user='taichi', processes=1)
    # # todo:
    # #   data_manager_jaxa_gsmap.get_raster_data(xxxx)
    #
    # # list_ftp_path = data_manager_jaxa_gsmap.get_ftp_path_list(
    # #     start_datetime=datetime.datetime(year=2020, month=1, day=26),
    # #     end_datetime=datetime.datetime(year=2020, month=4, day=26))
    # # print(list_ftp_path)
    # path_local = data_manager_jaxa_gsmap.download_from_ftp(
    #     path_ftp='/standard/GSMaP/3.GSMAP.M/03F/2016/GPMMRG_MAP_1601_M_L3S_MCM_03F.h5')
    # list_path_converted = data_manager_jaxa_gsmap.convert_from_hdf_to_gtiff(path_local,
    #                                                                         list_band_filter=[
    #                                                                             '//Grid/monthlyPrecipRate'])
    # path_converted = list_path_converted[0]
    # shapes = shapes_from_bbox(x_min=120.61, y_min=22.29, x_max=151.35, y_max=46.8, epsg_code=4326)
    # pick_part_raster(path_converted, os.path.join('data', os.path.basename(path_converted)), shapes)
    # # print(path_local, path_converted)

