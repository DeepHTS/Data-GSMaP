import os
import datetime
from tqdm import tqdm
import pandas as pd

import boto3

DICT_GSMAP_UNIT = {
    'Latitude': 'degrees',
    'Longitude': 'degrees',
    'hourlyPrecipRate': 'mm/hr',
    'hourlyPrecipRateGC': 'mm/hr',
    'monthlyPrecipRate': 'mm/hr',
    'monthlyPrecipRateGC': 'mm/hr',
    'standardDeviation': 'mm/hr',
    'gaugeQualityInfo': 'counts/day',
    'snowProbability': '%'
}


def get_meta_list_GSMaP(s3_bucket_name, s3_dir_parent):
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(s3_bucket_name)
    s3_cl = boto3.client('s3')
    bucket_location = s3_cl.get_bucket_location(Bucket=s3_bucket_name)

    list_path = []
    for bucket_object in bucket.objects.filter(Prefix=s3_dir_parent):
        list_path.append(bucket_object.key)

    list_dict_meta = []
    for path in tqdm(list_path):
        filename = os.path.basename(path)
        filename_part = filename.split('_')
        dict_info = {
            'filename': filename,
            'mission_id': filename_part[0],
            'sensor_id': filename_part[1],
            'scene_start_datetime': filename_part[2],
            'process_unit': filename_part[3],
            'process_level': filename_part[4],
            'product_version': filename_part[6],
            'bucket_name': s3_bucket_name,
            'path': path,
            'url': "https://{0}.s3-{1}.amazonaws.com/{2}".format(s3_bucket_name,
                                                                 bucket_location['LocationConstraint'], path)
        }
        if dict_info['process_unit'] == 'H':
            dict_info['scene_start_datetime'] = datetime.datetime.strptime(dict_info['scene_start_datetime'],
                                                                           '%y%m%d%H%M')
            dict_info['process_unit'] = 'Hourly'
        elif dict_info['process_unit'] == 'M':
            dict_info['scene_start_datetime'] = datetime.datetime.strptime(dict_info['scene_start_datetime'], '%y%m')
            dict_info['process_unit'] = 'Monthly'

        if dict_info['process_level'][-1] == 'S':
            dict_info['process_level'] = 'Standard'
        elif dict_info['process_level'][-1] == 'N':
            dict_info['process_level'] = 'Near real time'

        if len(filename_part) >= 9:
            dict_info['group_name'] = filename_part[7]
            dict_info['element'] = filename_part[8]
            dict_info['unit'] = DICT_GSMAP_UNIT.get(dict_info['element'], '')
        if len(filename_part) >= 11:
            dict_info['epsg_code'] = int(filename_part[9])
            dict_info['region_name'] = filename_part[10]

        list_dict_meta.append(dict_info)

    df = pd.DataFrame.from_dict(list_dict_meta)
    return df
