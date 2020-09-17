import os

from src.s3_metainfo import get_meta_list_GSMaP
from src.helper import transfer_to_s3

S3_BUCKET_NAME = 'xxxx'
FILE_NAME_HEAD = 'xxxx'


def main():
    s3_dir_parent = 'data/GSMaP/picked/standard/hourly'

    if not os.path.exists(s3_dir_parent):
        os.makedirs(s3_dir_parent)

    df = get_meta_list_GSMaP(S3_BUCKET_NAME, s3_dir_parent)
    path_local_csv = os.path.join(s3_dir_parent, FILE_NAME_HEAD + '.csv')
    df.to_csv(path_local_csv, index=False)
    print(path_local_csv)
    url_s3 = transfer_to_s3(path_local=path_local_csv,
                            dir_local_parent=s3_dir_parent,
                            dir_s3_parent=s3_dir_parent,
                            remove_local_file=False,
                            multiprocessing=False,
                            s3_bucket_name=S3_BUCKET_NAME)
    print(url_s3)


if __name__ == "__main__":
    main()
