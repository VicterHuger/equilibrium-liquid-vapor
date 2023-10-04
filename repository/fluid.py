import asyncio
from database.tortoise.models import Fluid
from database.excell.fluid_data import create_fluid_data


async def create_fluids():
    fluids_data = [create_fluid_data('Fluid1')]
    fluids_data.append(create_fluid_data('Fluid2'))
    fluids = await asyncio.gather(*[Fluid.create(**fluid) for fluid in fluids_data])
    for fluid in fluids:
        await fluid.save()
        print(f"Fluid {fluid.id}:")
        for field in Fluid._meta.fields_map:
            print(f"  {field}: {getattr(fluid, field)}")


async def delete_all_fluids():
    await Fluid.filter().delete()


async def get_all_fluids():
    return await Fluid.all()
