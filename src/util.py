import logging
import os
from pathlib import Path
from typing import Union

PathLike = Union[Path, str]

def setup_logging(logname : str, savepath : PathLike = os.getcwd()):
    logging.basicConfig(
        level=logging.INFO,  # Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(savepath, logname)), 
            logging.StreamHandler()
        ]
    )
    logging.info("Logger initialized")
