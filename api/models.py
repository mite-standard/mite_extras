from typing import Dict, List, Optional, Union

from pydantic import BaseModel


class ValidationItem(BaseModel):
    reaction_smarts: str
    expected_positive: List[tuple[str,str]]
    expected_negative: List[tuple[str,str]]


class Item(BaseModel):
    validation: ValidationItem = ValidationItem()
    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "validation": {
                        "reaction_smarts": "foo>>bar",
                        "expected_positive": [["foo", "bar"]],
                        "expected_negative": [["alice", "bob"]],
                        "option": {
                            "date_min": "1999-12-31",
                            "date_max": "2024-01-01",
                            "journal": "Natural Product Research",
                        },
                    },
                }
            ]
        }
    }


class ValidationObject(BaseModel):
    result: bool


class ValidationResult(BaseModel):
    ids: List[int]
    objects: Optional[Dict[int, ValidationObject]] = None
