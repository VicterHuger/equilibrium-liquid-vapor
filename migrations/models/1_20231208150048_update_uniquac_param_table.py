from tortoise import BaseDBAsyncClient


async def upgrade(db: BaseDBAsyncClient) -> str:
    return """
        ALTER TABLE "uniquacparm" RENAME COLUMN "param" TO "param_0";
        ALTER TABLE "uniquacparm" ADD COLUMN "param_t" DOUBLE PRECISION NOT NULL  DEFAULT 0;
        """


async def downgrade(db: BaseDBAsyncClient) -> str:
    return """
        ALTER TABLE "uniquacparm" RENAME COLUMN "param_0" TO "param";
        ALTER TABLE "uniquacparm" DROP COLUMN "param_t";
        """
