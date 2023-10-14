"""Repository file for peng robinson params creating, deleting and getting operations
"""
from database.tortoise.models import PengRobinsonParam
from database.seed_data.peng_robinson_params_data import PENG_ROBINSON_PARAMS_VALUES
from repository.component import get_component


async def create_peng_robinson_params():
    """A function to create peng robinson params in datbase
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
    """A function to delete all peng robinson params in datbase
    """
    await PengRobinsonParam.filter().delete()


async def get_all_peng_robinson_params():
    """A function to get all peng robinson params in datbase
    """
    return await PengRobinsonParam.all()
