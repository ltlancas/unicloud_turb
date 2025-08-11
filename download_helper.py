## download_helper.py
## Code to help download the default velocity cube from where it is stored on zenodo
## The normal storage location is https://zenodo.org/records/16795651


import requests
from tqdm import tqdm

def downloader(fname, url, message):
    """
    Method for downloading data from the web.

    Args:
        fname (str or Path): The file name or path to save the downloaded file.
        url (str): The URL of the file to download.
        message (str): Optional message to display before downloading.

    Raises:
        Exception: If the download fails.

    """
    if message is not None:
        print(message)

    try:
        response = requests.get(url, stream=True, timeout=10)
    except requests.exceptions.Timeout:
        raise Exception('Request timed out. Check your internet connection.')
    except requests.exceptions.RequestException as e:
        raise Exception(f'Request failed: {e}')

    # Get file size
    total_size = int(response.headers.get('content-length', 0))
    # create a progress bar
    tqdm_args = dict(desc='Downloading', total=total_size, unit='B',
                    unit_scale=True, unit_divisor=1024)
    # write the file
    with open(fname, 'wb') as f, tqdm(**tqdm_args) as prog_bar:
        for chunk in response.iter_content(chunk_size=1024):
            f.write(chunk)
            prog_bar.update(len(chunk))
