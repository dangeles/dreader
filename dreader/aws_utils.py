import boto3
import os
import sys

def download_from_s3(bucket_name, prefix, download_path):
    """
    A function to download files from an S3 bucket.
    
    Parameters
    ----------
    bucket_name: str
        Name of the S3 bucket.
    prefix: str
        Prefix to filter files in the bucket.
    download_path: str
        Path where files will be downloaded
    """
    s3 = boto3.client('s3')
    paginator = s3.get_paginator('list_objects_v2')
    for page in paginator.paginate(Bucket=bucket_name, Prefix=prefix):
        if 'Contents' in page:
            for obj in page['Contents']:
                key = obj['Key']
                local_path = os.path.join(download_path, key[len(prefix):])
                local_dir = os.path.dirname(local_path)
                if not os.path.exists(local_dir):
                    os.makedirs(local_dir)
                if not key.endswith('/'):
                    s3.download_file(bucket_name, key, local_path)

def fetch_data(dataset, sample=None, directory=None):
    """
    A function to download desired datasets from organogenesis-data-processed S3 bucket.
    
    Parameters
    ----------
    dataset: str
        Name of the dataset.
    sample: str
        Name of the sample to download. If None, the entire dataset will be downloaded.
    directory: str
        Path to download data to. If None, downloads to current working directory
    """
    bucket_name = 'organogenesis-data-processed'
    if sample:
        prefix = f'{dataset}/{sample}/'
    else:
        prefix = f'{dataset}/'
    
    if directory == None:
        download_path = os.path.join(os.getcwd(), prefix)
    else:
        download_path = os.path.join(directory, prefix)
    
    print('Downloading data...')    
    download_from_s3(bucket_name, prefix, download_path)
    print(f'Download completed: {download_path}')
