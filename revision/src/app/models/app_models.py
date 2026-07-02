from typing import Optional

from pydantic import BaseModel, ConfigDict


# Define UserModel
class UserModel(BaseModel):
    # Allow model extension
    model_config = ConfigDict(extra="allow")

    username: str
    email: Optional[str] = None
