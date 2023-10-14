"""file to inialize project"""
from dotenv import load_dotenv
from tortoise import Tortoise, run_async
from database.tortoise.tortoise_config import TORTOISE
from seed import seed

load_dotenv()


async def create_tables():
    """function to create tortoise tables in database
    """
    await Tortoise.init(
        config=TORTOISE,
        modules={"models": ["database.tortoise.models"]}
    )
    await Tortoise.generate_schemas()


async def init():
    """function to inialize project
    """
    await create_tables()
    print('Tables Created')

    is_seed_requested = input(
        "Do you want to populate database with seed data? Answer with 1 for 'YES' or 2 for 'NO'\n")
    if is_seed_requested == '1':
        await seed()
    await Tortoise.close_connections()

run_async(init())
