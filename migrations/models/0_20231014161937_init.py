from tortoise import BaseDBAsyncClient


async def upgrade(db: BaseDBAsyncClient) -> str:
    return """
        CREATE TABLE IF NOT EXISTS "component" (
    "id" SERIAL NOT NULL PRIMARY KEY,
    "name" VARCHAR(20) NOT NULL UNIQUE,
    "critical_temperature" DOUBLE PRECISION NOT NULL,
    "critical_pressure" DOUBLE PRECISION NOT NULL,
    "acentric_factor" DOUBLE PRECISION NOT NULL,
    "aat1" DOUBLE PRECISION NOT NULL  DEFAULT 1,
    "aat2" DOUBLE PRECISION NOT NULL  DEFAULT 1,
    "aat3" DOUBLE PRECISION NOT NULL  DEFAULT 1,
    "r" DOUBLE PRECISION NOT NULL  DEFAULT 1,
    "q" DOUBLE PRECISION NOT NULL  DEFAULT 1,
    "created_at" TIMESTAMPTZ NOT NULL  DEFAULT CURRENT_TIMESTAMP,
    "updated_at" TIMESTAMPTZ NOT NULL  DEFAULT CURRENT_TIMESTAMP
);
COMMENT ON TABLE "component" IS 'The Component Model';
CREATE TABLE IF NOT EXISTS "fluid" (
    "id" SERIAL NOT NULL PRIMARY KEY,
    "name" VARCHAR(50) NOT NULL UNIQUE,
    "temp_pressure_coordinates" JSONB,
    "created_at" TIMESTAMPTZ NOT NULL  DEFAULT CURRENT_TIMESTAMP,
    "updated_at" TIMESTAMPTZ NOT NULL  DEFAULT CURRENT_TIMESTAMP
);
COMMENT ON COLUMN "fluid"."temp_pressure_coordinates" IS 'Should be added an array of dictionaries with         temperature (K) and pressure (bar) properties';
COMMENT ON TABLE "fluid" IS 'The Fluid Model';
CREATE TABLE IF NOT EXISTS "componentfluid" (
    "id" SERIAL NOT NULL PRIMARY KEY,
    "composition" DOUBLE PRECISION NOT NULL  DEFAULT 0,
    "component_id" INT NOT NULL REFERENCES "component" ("id") ON DELETE CASCADE,
    "fluid_id" INT NOT NULL REFERENCES "fluid" ("id") ON DELETE CASCADE,
    CONSTRAINT "uid_componentfl_compone_1fd3df" UNIQUE ("component_id", "fluid_id")
);
COMMENT ON TABLE "componentfluid" IS 'ComponentFluid represents the many-to-many relationship with an additional composition field.';
CREATE TABLE IF NOT EXISTS "pengrobinsonparam" (
    "id" SERIAL NOT NULL PRIMARY KEY,
    "param" DOUBLE PRECISION NOT NULL  DEFAULT 0,
    "componentA_id" INT NOT NULL REFERENCES "component" ("id") ON DELETE CASCADE,
    "componentB_id" INT NOT NULL REFERENCES "component" ("id") ON DELETE CASCADE,
    CONSTRAINT "uid_pengrobinso_compone_ed2ca9" UNIQUE ("componentA_id", "componentB_id")
);
COMMENT ON TABLE "pengrobinsonparam" IS 'PengRobinsonParam for a combination of two components';
CREATE TABLE IF NOT EXISTS "uniquacparm" (
    "id" SERIAL NOT NULL PRIMARY KEY,
    "param" DOUBLE PRECISION NOT NULL  DEFAULT 0,
    "componentA_id" INT NOT NULL REFERENCES "component" ("id") ON DELETE CASCADE,
    "componentB_id" INT NOT NULL REFERENCES "component" ("id") ON DELETE CASCADE,
    CONSTRAINT "uid_uniquacparm_compone_3a2da6" UNIQUE ("componentA_id", "componentB_id")
);
COMMENT ON TABLE "uniquacparm" IS 'Uniquac Param model for two components';
CREATE TABLE IF NOT EXISTS "aerich" (
    "id" SERIAL NOT NULL PRIMARY KEY,
    "version" VARCHAR(255) NOT NULL,
    "app" VARCHAR(100) NOT NULL,
    "content" JSONB NOT NULL
);"""


async def downgrade(db: BaseDBAsyncClient) -> str:
    return """
        """
