""" Tortoise configuration with postgres sql   """
import os
from dotenv import load_dotenv

load_dotenv()

TORTOISE = {
    "connections": {
        "default": os.getenv("DATABASE_URL"),
    },
    "apps": {
        "models": {
            "models": ["database.tortoise.models", "aerich.models"],
            "default_connection": "default",
        }
    },
    "use_tz": True,
    "timezone": "UTC",
    "generate_schemas": True
}
