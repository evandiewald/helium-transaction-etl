from sqlalchemy.orm import declarative_base, relationship
from sqlalchemy import Column, Text, BigInteger, Integer, Boolean, Float, DateTime, ForeignKey, CheckConstraint, Enum
from sqlalchemy.dialects.postgresql import UUID, JSONB
from geoalchemy2 import Geometry
import uuid
import os
import enum


Base = declarative_base(bind=os.getenv("POSTGRES_CONNECTION_STR"))


class witness_invalid_reason_type(enum.Enum):
    witness_rssi_too_high = 1
    incorrect_frequency = 2
    witness_not_same_region = 3
    witness_too_close = 4
    witness_on_incorrect_channel = 5
    witness_too_far = 6


class payment_type(enum.Enum):
    payment_v1 = 1
    payment_v2 = 2


class ChallengeReceiptsParsed(Base):
    __tablename__ = "challenge_receipts_parsed"

    block = Column(BigInteger, nullable=False, index=True)
    hash = Column(Text, nullable=False, primary_key=True)
    time = Column(BigInteger, nullable=False)
    challenger = Column(Text, ForeignKey("gateway_inventory.address"), nullable=False)
    transmitter_address = Column(Text, ForeignKey("gateway_inventory.address"), nullable=False)
    tx_power = Column(Integer)
    origin = Column(Text)
    witness_address = Column(Text, ForeignKey("gateway_inventory.address"), primary_key=True)
    witness_is_valid = Column(Boolean, index=True)
    witness_invalid_reason = Column(Enum(witness_invalid_reason_type))
    witness_signal = Column(Integer)
    witness_snr = Column(Float)
    witness_channel = Column(Integer)
    witness_datarate = Column(Text)
    witness_frequency = Column(Float)
    witness_timestamp = Column(BigInteger)
    ra = Column(Float)
    rq = Column(Float)
    rv = Column(Float)
    rz = Column(Float)
    rsk = Column(Float)
    rku = Column(Float)
    deepest_barrier = Column(Float)
    n_barriers = Column(Integer)
    distance_km = Column(Float, index=True)


class PaymentsParsed(Base):
    __tablename__ = "payments_parsed"

    block = Column(BigInteger)
    hash = Column(Text, primary_key=True)
    time = Column(BigInteger)
    payer = Column(Text, index=True)
    payee = Column(Text, primary_key=True)
    amount = Column(BigInteger)
    type = Column(Enum(payment_type))
    fee = Column(BigInteger)
    nonce = Column(BigInteger)


class DataCredits(Base):
    __tablename__ = "data_credits"

    block = Column(BigInteger)
    hash = Column(Text, primary_key=True)
    client = Column(Text, ForeignKey("gateway_inventory.address"), primary_key=True)
    num_dcs = Column(Integer)
    num_packets = Column(Integer)


class GatewayInventory(Base):
    __tablename__ = "gateway_inventory"

    address = Column(Text, primary_key=True, nullable=False, unique=True)
    owner = Column(Text)
    location = Column(Text)
    last_poc_challenge = Column(BigInteger)
    last_poc_onion_key_hash = Column(Text)
    first_block = Column(BigInteger)
    last_block = Column(BigInteger)
    nonce = Column(BigInteger)
    name = Column(Text)
    first_timestamp = Column(DateTime)
    reward_scale = Column(Float)
    elevation = Column(Integer)
    gain = Column(Integer)
    location_hex = Column(Text)
    mode = Column(Text)
    payer = Column(Text, index=True)


class Denylist(Base):
    __tablename__ = "denylist"

    address = Column(Text, ForeignKey("gateway_inventory.address"), primary_key=True)


class FollowerInfo(Base):
    __tablename__ = "follower_info"

    name = Column(Text, primary_key=True)
    value = Column(BigInteger)
