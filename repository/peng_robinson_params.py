"""Repository file for peng robinson params creating, deleting and getting operations
"""
from database.tortoise.models import PengRobinsonParam
from database.seed_data.peng_robinson_params_data import PENG_ROBINSON_PARAMS_VALUES
from repository.component import get_component
from tortoise.exceptions import DoesNotExist


async def create_peng_robinson_params():
    """A function to create peng robinson params in database
    """

    peng_robinson_params_data = PENG_ROBINSON_PARAMS_VALUES.items()
    for component_a_name, param_dict in peng_robinson_params_data:
        component_a_obj = await get_component(component_a_name)
        for component_b_name, param_value in param_dict.items():
            component_b_obj = await get_component(component_b_name)
            peng_robinson_param = PengRobinsonParam(
                componentA=component_a_obj, componentB=component_b_obj, param=param_value)
            await peng_robinson_param.save()


async def delete_all_peng_robinson_params():
    """A function to delete all peng robinson params in database
    """
    await PengRobinsonParam.filter().delete()


async def get_all_peng_robinson_params():
    """A function to get all peng robinson params in database
    """
    return await PengRobinsonParam.all()


async def get_peng_robinson_param_by_component_id(component_a_id: str, component_b_id: str):
    """A function to get peng robinson binary param for two components ids provided
    """
    try:
        return await PengRobinsonParam.get(componentA_id=component_a_id, componentB_id=component_b_id)
    except DoesNotExist:
        try:
            return await PengRobinsonParam.get(componentB_id=component_a_id, componentA_id=component_b_id)
        except Exception as err:
            raise Exception(err)
