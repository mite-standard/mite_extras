import logging
from datetime import datetime
from typing import Any

from fastapi import HTTPException, status

from api.models import (
    Item,
    ValidationItem,
)
from model.data_model import DataModel

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def get_validation_for_reaction_smarts(item: Item, dm: DataModel) -> dict[int, str]:
    # TODO
    # return dm.todo(smarts)
