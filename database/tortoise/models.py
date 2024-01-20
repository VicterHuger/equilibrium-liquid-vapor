"""
Tortoise models   
"""
from tortoise import fields
from tortoise.models import Model

MODEL_COMPONENT = "models.Component"


class Component(Model):
    """
    The Component Model

    """
    id = fields.IntField(pk=True)
    name = fields.CharField(max_length=20, unique=True, null=False)
    critical_temperature = fields.FloatField(null=False)
    critical_pressure = fields.FloatField(null=False)
    acentric_factor = fields.FloatField(null=False)
    aat1 = fields.FloatField(default=1., null=False)
    aat2 = fields.FloatField(default=1., null=False)
    aat3 = fields.FloatField(default=1., null=False)
    r = fields.FloatField(default=1.,)
    q = fields.FloatField(default=1.,)
    created_at = fields.DatetimeField(auto_now_add=True)
    updated_at = fields.DatetimeField(auto_now=True)

    def __str__(self):
        return self.name


class Fluid(Model):
    """
    The Fluid Model

    """
    id = fields.IntField(pk=True)
    name = fields.CharField(max_length=50, unique=True, null=False)
    temp_pressure_coordinates = fields.JSONField(
        null=True, description="Should be added an array of dictionaries with \
        temperature (K) and pressure (bar) properties")
    created_at = fields.DatetimeField(auto_now_add=True)
    updated_at = fields.DatetimeField(auto_now=True)


class ComponentFluid(Model):
    """
    ComponentFluid represents the many-to-many relationship with an additional composition field.

    """

    id = fields.IntField(pk=True)
    component = fields.ForeignKeyField(
        MODEL_COMPONENT, related_name="component_fluids")
    fluid = fields.ForeignKeyField(
        "models.Fluid", related_name="fluid_components")
    composition = fields.FloatField(default=0., null=False)

    class Meta:
        """Meta class with unique_together is used to ensure 
            that the combination of component and fluid  
            is unique in the ComponentFluid table.
        """
        unique_together = ("component", "fluid")


class PengRobinsonParam(Model):
    """
    PengRobinsonParam for a combination of two components

    """

    id = fields.IntField(pk=True)
    componentA = fields.ForeignKeyField(
        MODEL_COMPONENT, related_name="peng_robinson_param_A")
    componentB = fields.ForeignKeyField(
        MODEL_COMPONENT, related_name="peng_robinson_param_B")
    param = fields.FloatField(default=0., null=False)

    class Meta:
        """Meta class with unique_together is used to ensure 
            that the combination of componentA and componentB  
            is unique in the PengRobinsonParam table.
        """
        unique_together = ("componentA", "componentB")


class UniquacParm(Model):
    """
    Uniquac Params model for two components: 
    Uij_0: List[List[float]]
    Uij_T: [List[float]]
    """
    id = fields.IntField(pk=True)
    componentA = fields.ForeignKeyField(MODEL_COMPONENT, 'uniquac_param_A')
    componentB = fields.ForeignKeyField(MODEL_COMPONENT, 'uniquac_param_B')
    param_0 = fields.FloatField(null=False, default=0.0)
    param_t = fields.FloatField(null=False, default=0.0)

    class Meta:
        """Meta class with unique_together is used to ensure 
        that the combination of componentA and componentB  
        is unique in the UniquacParam table.
        """
        unique_together = ("componentA", "componentB")

    async def save(self, *args, **kwargs):
        if self.componentA == self.componentB:
            self.param_0 = 0.0
            self.param_t = 0.0
        super().save(*args, **kwargs)
