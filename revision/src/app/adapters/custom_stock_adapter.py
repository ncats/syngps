import os
import pandas as pd
from syngps.adapters import AbstractInventoryAdapter
from app.logging_config import logger
from syngps.models.models import Availability


class CustomStockInventoryAdapter(AbstractInventoryAdapter):
    def __init__(self, stock_inv_file_path, inchikey_col="inchikey"):
        # Inv file path is required
        if not stock_inv_file_path:
            raise ValueError("Stock inventory file path is required.")
        
        self.inchikey_col = inchikey_col
        if not os.path.exists(stock_inv_file_path):
            logger.warning(f"Stock inventory file not found at {stock_inv_file_path}. Loading empty inventory.")
            self.inventory_df = pd.DataFrame(columns=[inchikey_col])
        else:
            self.inventory_df = pd.read_csv(stock_inv_file_path, header=None, names=[inchikey_col])
        self.inchikey_set = set(self.inventory_df[inchikey_col].dropna().unique())
        self.count = len(self.inchikey_set)

        logger.info(f"CustomStockInventoryAdapter initialized with {self.count} unique inchikeys from {stock_inv_file_path}")

    def in_inventory_by_inchikey(self, inchikey: str) -> bool:
        return inchikey in self.inchikey_set
    
    def inchikey_inventory_status(self, inchikeys: list) -> dict:
        return {inchikey: self.in_inventory_by_inchikey(inchikey) for inchikey in inchikeys}
    
    def adapter_status(self) -> dict:
        return {
            "custom_stock_inventory_adapter": {
                "status": "connected" if self.inventory_df is not None else "not connected",
                "num_compounds": self.count if self.inventory_df is not None else 0,
            }
        }
    
    def inventory_status_details(self, inchikey: str) -> Availability:
        available = self.in_inventory_by_inchikey(inchikey)
        return Availability(
            inchikey=inchikey,
            inventory={"available": available},
            commercial_availability=None
        )