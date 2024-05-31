import logging
from contextlib import asynccontextmanager

from fastapi import Depends, FastAPI
from fastapi_versioning import VersionedFastAPI, version

from api.models import (
    Item,
    ValidationObject,
    ValidationResult,
)
# from api.queries import (
#     todo,
# )
from model.data_model import DataModel

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

description = """
MITEFast API helps you do awesome stuff. 🚀
"""

data_model: None | DataModel = None


def get_data_model() -> DataModel:
    """
    A bit messy, but that way we can inject our own in tests.
    (See lotus-search)
    :return:
    """
    global data_model
    return data_model


app = FastAPI(
    title="MITE FastAPI",
    description=description,
    summary="An awesome way to validate reaction SMARTS.",
    version="1.0",
    # TODO
    # terms_of_service="http://example.com/terms/",
    # contact={
    #     "name": "Deadpoolio the Amazing",
    #     "url": "http://x-force.example.com/contact/",
    #     "email": "dp@x-force.example.com",
    # },
    license_info={
        "name": "Apache 2.0",
        "url": "https://www.apache.org/licenses/LICENSE-2.0.html",
    },
)


@asynccontextmanager
async def lifespan(_: FastAPI):
    global data_model
    data_model = DataModel()
    yield
    data_model = None


@app.post("/validate")
@version(1, 0)
async def todo(
    item: Item, dm: DataModel = Depends(get_data_model)
) -> ValidationResult:
    # TODO


LOGGING_CONFIG = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "default": {
            "()": "logging.Formatter",
            "fmt": "%(levelname)s %(name)s@%(lineno)d %(message)s",
        },
    },
    "handlers": {
        "default": {
            "formatter": "default",
            "class": "my_project.ColorStreamHandler",
            "stream": "ext://sys.stderr",
        },
    },
    "loggers": {
        "": {"handlers": ["default"], "level": "TRACE"},
    },
}

app = VersionedFastAPI(
    app, enable_latest=True, log_config=LOGGING_CONFIG, lifespan=lifespan
)
