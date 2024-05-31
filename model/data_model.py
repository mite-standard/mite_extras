import logging

from api.models import ValidationObject

logging.basicConfig()
log = logging.getLogger()
log.setLevel(logging.WARNING)
requests_log = logging.getLogger("requests.packages.urllib3")
requests_log.setLevel(logging.WARNING)
requests_log.propagate = True


class DataModel:
    # TODO

    # def todo(self, smarts: int) -> dict | None:
    #     return self.todo(smarts)
