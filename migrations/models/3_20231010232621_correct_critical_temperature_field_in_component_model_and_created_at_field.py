from tortoise import BaseDBAsyncClient


async def upgrade(db: BaseDBAsyncClient) -> str:
    return """
        ALTER TABLE "component" RENAME COLUMN "critical_temperatue" TO "critical_temperature";
        ALTER TABLE "component" ALTER COLUMN "created_at" TYPE TIMESTAMPTZ USING "created_at"::TIMESTAMPTZ;"""


async def downgrade(db: BaseDBAsyncClient) -> str:
    return """
        ALTER TABLE "component" RENAME COLUMN "critical_temperature" TO "critical_temperatue";
        ALTER TABLE "component" ALTER COLUMN "created_at" TYPE TIMESTAMPTZ USING "created_at"::TIMESTAMPTZ;"""
