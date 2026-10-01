import pandas as pd
import numpy as np
import os
import pathlib
from typing import List, Optional
import pyarrow as pa
import pyarrow.parquet as pq

def dump_pandas(
    df: pd.DataFrame, 
    fname: str,
    location: str, 
    subdirs: Optional[List[str]] = None,
    pyarrow: bool = False
) -> None:
    """
    Dump a pandas dataframe to disk at location/'/'.join(subdirs)/fname.
    """

    # Set destination
    subdirs = subdirs or []
    merged_subdirs = os.path.sep.join(subdirs)
    destination = (os.path.join(location, os.path.join(merged_subdirs, fname)))

    # Check existence of destination and create if necessary
    dirname = os.path.dirname(destination)
    if not os.path.exists(dirname):
        pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)

    # Save to parquet and check
    if os.path.exists(destination):
        os.remove(destination)

    if not pyarrow:
        df.to_parquet(destination)
    else:
        pq.write_table(pa.Table.from_pandas(df), destination)
    assert os.path.isfile(destination), "Saved parquet does not exist."
