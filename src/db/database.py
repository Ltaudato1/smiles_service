from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base
from dotenv import load_dotenv
from os import getenv

load_dotenv()

user = getenv('POSTGRES_USER')
password = getenv('POSTGRES_PASSWORD')
db = getenv('POSTGRES_DB')

DATABASE_URL = f"postgresql://{user}:{password}@db:5432/{db}"

engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()
