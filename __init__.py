"""file to inialize project"""
from dotenv import load_dotenv
from tortoise import Tortoise, run_async, exceptions
from database.tortoise.tortoise_config import TORTOISE
from service import bubble_point, dew_point
from seed import seed
from repository import fluid, composition, peng_robinson_params
from utils import binary_model_param, mixing_rule_param, uniquac_properties, component_properties, alpha_model
from typing import List
import math

load_dotenv()


async def create_tables():
    """function to create tortoise tables in database
    """
    await Tortoise.init(
        config=TORTOISE,
        modules={"models": ["database.tortoise.models"]}
    )
    await Tortoise.generate_schemas()

# def calculate_pseudo_component_properties(cricital_temperature:float, critical_pressure:float, acentric_factor:float, specific_gravity:float, molecular_weight:float):
#     ck =


async def init():
    """function to inialize project
    """

    try:
        await create_tables()
        print('Tables Created \n')

        is_seed_requested = input(
            "Do you want to populate database with seed data? Answer with 1 for 'YES' or 2 for 'NO'\n")
        if is_seed_requested == '1':
            await seed()

        binary_model_params_input = input(
            '\n Which binary model param would you like to use?  Answer 1 for "Peng Robinson" or 2 for "Peng Robinson + UNIQUAC"\n')
        binary_model_params_choice = binary_model_param.BinaryModelParam.PR_UNIQUAC if binary_model_params_input == '2' else binary_model_param.BinaryModelParam.PR

        alpha_model_input = input(
            '\n Which alpha model function would you like to use?  Answer 1 for "Peng Robinson - PR78" or 2 for "Almeida-Aznar- Telles (AAT) function"\n')
        alpha_model_choice = alpha_model.AlphaModel.AAT if alpha_model_input == '2' else alpha_model.AlphaModel.PR87

        alpha_aat_params: List[List[float]] | None = None

        mixing_rule_params_choice: mixing_rule_param.MixingRuleParam | None = None
        uniquac_properties_data: uniquac_properties.UniquacProperties = uniquac_properties.UniquacProperties([
        ], [], [], [])

        if binary_model_params_choice == binary_model_param.BinaryModelParam.PR_UNIQUAC:
            mixing_rule_params_input = input(
                '\n Which mixing rule params would you like to use?  Answer 1 for "Heidemann Kokal - HK" or 2 for "LCVM"\n')
            mixing_rule_params_choice = mixing_rule_param.MixingRuleParam.LCVM if mixing_rule_params_input == '2' else mixing_rule_param.MixingRuleParam.HK
            # TODO: Implement repository for GETTING UNIQUAC PARAMS and UPDATE uniquac_properties_data

        fluid_input_choice_options_text = ", ".join([fluid_db.name for fluid_db in (await fluid.get_all_fluids())])
        fluid_choice = input(
            f"\n Choose a fluid:\n{fluid_input_choice_options_text} \n")

        try:
            fluid_db = await fluid.get_fluid(fluid_name=fluid_choice)
        except exceptions.DoesNotExist:
            raise exceptions.DoesNotExist(
                f'Fluid name provided "{fluid_choice}" is not a valid fluid name! Valid options: {fluid_input_choice_options_text}')

        components_fluids = await composition.get_components_fluids_by_fluid(fluid=fluid_db)

        compostions = [
            component_fluid.composition for component_fluid in components_fluids]

        acentric_factors: List[float] = []
        critical_temperatures: List[float] = []
        critical_pressures: List[float] = []
        peng_robinson_binary_params: List[List[float]] = []

        if alpha_model_choice == alpha_model.AlphaModel.AAT:
            alpha_aat_params = [[], [], []]
            for component_fluid in components_fluids:
                component = component_fluid.component.__dict__

                acentric_factors.append(component['acentric_factor'])
                critical_temperatures.append(component['critical_temperature'])
                critical_pressures.append(component['critical_pressure'])
                alpha_aat_params[0].append(component['aat1'])
                alpha_aat_params[1].append(component['aat2'])
                alpha_aat_params[2].append(component['aat3'])

                peng_robinson_binary_params.append([])
                peng_robinson_binary_params_length = peng_robinson_binary_params.__len__()

                for other_component_fluid in components_fluids:
                    other_component = other_component_fluid.component.__dict__
                    peng_robinson_param = await peng_robinson_params.get_peng_robinson_param_by_component_id(component_a_id=component['id'], component_b_id=other_component['id'])
                    peng_robinson_binary_params[peng_robinson_binary_params_length-1].append(
                        peng_robinson_param.param)

        else:
            for component_fluid in (components_fluids):
                component = component_fluid.component.__dict__

                acentric_factors.append(component['acentric_factor'])
                critical_temperatures.append(component['critical_temperature'])
                critical_pressures.append(component['critical_pressure'])

                peng_robinson_binary_params.append([])
                peng_robinson_binary_params_length = peng_robinson_binary_params.__len__()

                for other_component_fluid in components_fluids:
                    other_component = other_component_fluid.component.__dict__
                    peng_robinson_param = await peng_robinson_params.get_peng_robinson_param_by_component_id(component_a_id=component['id'], component_b_id=other_component['id'])
                    peng_robinson_binary_params[peng_robinson_binary_params_length-1].append(
                        peng_robinson_param.param)

        critical_acentric_properties = component_properties.CriticalAcentricProperties(
            acentric_factors, critical_temperatures, critical_pressures)
        aat_params = alpha_model.AATParams(
            aat1=alpha_aat_params[0], aat2=alpha_aat_params[1], aat3=alpha_aat_params[2]) if alpha_model_choice == alpha_model.AlphaModel.AAT else None
        component_properties_data = component_properties.ComponentProperties(
            peng_robinson_binary_params=peng_robinson_binary_params, critical_acentric_properties=critical_acentric_properties, alpha_model=alpha_model_choice, aat_params=aat_params)

        s = [[str(e) for e in row]
             for row in peng_robinson_binary_params]
        lens = [max(map(len, col)) for col in zip(*s)]
        fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
        table = [fmt.format(*row) for row in s]
        print('Peng Robinson Params')
        print('\n'.join(table))

        # for i in range(fluid_db.temp_pressure_coordinates.__len__()):
        #     bubble_pressure = bubble_point.calc_bubble_pressure(
        #         number_components=components_fluids.__len__(),
        #         liquid_compositions=compostions,
        #         temperature=fluid_db.temp_pressure_coordinates[i]['temperature'],
        #         binary_model_param=binary_model_params_choice,
        #         component_properties=component_properties_data,
        #         uniquac_properties=uniquac_properties_data,
        #         mixing_rule_param=mixing_rule_params_choice,
        #     )
        #     print('buble_pressure', bubble_pressure)

        for i in range(fluid_db.temp_pressure_coordinates.__len__()):
            dew_pressure = dew_point.calc_dew_pressure(
                number_components=components_fluids.__len__(),
                vapor_compositions=compostions,
                temperature=fluid_db.temp_pressure_coordinates[i]['temperature'],
                binary_model_param=binary_model_params_choice,
                component_properties=component_properties_data,
                uniquac_properties=uniquac_properties_data,
                mixing_rule_param=mixing_rule_params_choice,
            )

            print('\n Dew pressure expected',
                  fluid_db.temp_pressure_coordinates[i]['pressure'])
            print('\n Dew pressure calculated', dew_pressure)

            print('\n DABS Error (%):', abs(
                dew_pressure-fluid_db.temp_pressure_coordinates[i]['pressure'])/fluid_db.temp_pressure_coordinates[i]['pressure']*100)

    except KeyboardInterrupt:
        print('\n Calculations aborted by user')

    except Exception as err:
        print('Error occured in init function:', err)
        raise Exception(err)

    finally:
        await Tortoise.close_connections()

run_async(init())
