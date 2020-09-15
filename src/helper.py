from tqdm import tqdm
from multiprocessing import Pool
import boto3
import os

def argwrapper(args):
    return args[0](*args[1:])

def imap_unordered_bar(func, args, n_processes=15, extend=False, tqdm_disable=False,
                       init=None, credentials=None):
    """ execute for loop with multiprocessing

    Args:
        func (def): function of arg wrapper (argwrapper)
        args (list of args):  ex. [(target function, args of functions) for xx in xxxx]
        n_processes (int): number of processes
        django_process (bool): True for connection DB through django
        extend (bool): If True, each results will be merged by xxxx.extend(list)
        tqdm_disable (bool): If True, tqdm bar will not shown
        init (func): look def init(), call init(*credentials) before use this func
        credentials (list): [ftp_address, user name, password]

    Returns: appended or extended list

    """
    if (init is not None) and (credentials is not None):
        p = Pool(n_processes, initializer=init, initargs=credentials)
    else:
        p = Pool(n_processes)
    res_list = []
    with tqdm(total=len(args), disable=tqdm_disable) as pbar:
        for i, res in tqdm(enumerate(p.imap_unordered(func, args))):
            pbar.update()
            if extend:
                res_list.extend(res)
            else:
                res_list.append(res)
    pbar.close()
    p.close()
    p.join()
    return res_list

def transfer_to_s3(path_local, dir_local_parent=None, dir_s3_parent=None, remove_local_file=False, multiprocessing=False,
                   s3_bucket_name=None):
    """ transfer local file to s3 bucket

    Args:
        path_local (char): path to a target file on local (ex. 'macro_yield/data/converted/jaxa/GCOM-C/xxxxx' )
        dir_local_parent (char): parent directory for getting child path (ex. 'macro_yield/data/converted')
        dir_s3_parent (char): path on S3 bucket. The child path is attached after this path (ex. 'Macro_Yield/data/converted')
        remove_local_file (bool): If True, path_local will be removed
        multiprocessing (bool): If True, this code can be suit to multiprocessing

    Returns: url of the saved file on S3

    """
    if not multiprocessing:
        s3 = boto3.resource('s3')
        bucket = s3.Bucket(s3_bucket_name)
        s3 = boto3.client('s3')
        bucket_location = s3.get_bucket_location(Bucket=s3_bucket_name)
    else:
        session = boto3.session.Session()
        s3 = session.resource('s3')
        bucket = s3.Bucket(s3_bucket_name)
        s3 = session.client('s3')
        bucket_location = s3.get_bucket_location(Bucket=s3_bucket_name)

    if dir_local_parent is not None:
        path_child = path_local.split(dir_local_parent)[1][1:]
    else:
        path_child = str(path_local)
    if dir_s3_parent is not None:
        dest_path = os.path.join(dir_s3_parent, path_child)
    else:
        dest_path = str(path_child)
    bucket.upload_file(path_local, dest_path)
    url_s3 = "https://{0}.s3-{1}.amazonaws.com/{2}".format(
        s3_bucket_name,
        bucket_location['LocationConstraint'],
        dest_path)

    if remove_local_file:
        os.remove(path_local)

    return url_s3