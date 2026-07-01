import pandas as pd
from syngps.adapters import AbstractInventoryAdapter
from app.logging_config import logger
from syngps.models.models import Availability


class AskcosInventoryAdapter(AbstractInventoryAdapter):
    def __init__(self, askcos_inv_file_path, inchikey_col="inchikey"):
        # Inv file path is required
        if not askcos_inv_file_path:
            raise ValueError("Askcos inventory file path is required.")
        
        self.inchikey_col = inchikey_col
        # Split based off of tsv vs csv
        if askcos_inv_file_path.endswith(".tsv"):
            self.inventory_df = pd.read_csv(askcos_inv_file_path, sep="\t")
        else:
            self.inventory_df = pd.read_csv(askcos_inv_file_path)
        self.inchikey_set = set(self.inventory_df[inchikey_col].dropna().unique())
        self.count = len(self.inchikey_set)
        self.source_column = "source" if "source" in self.inventory_df.columns else None
        self.url_column = "URL" if "URL" in self.inventory_df.columns else None

        logger.info(f"AskcosInventoryAdapter initialized with {self.count} unique inchikeys from {askcos_inv_file_path}")

    def in_inventory_by_inchikey(self, inchikey: str) -> bool:
        return inchikey in self.inchikey_set
    
    def inchikey_inventory_status(self, inchikeys: list) -> dict:
        return {inchikey: self.in_inventory_by_inchikey(inchikey) for inchikey in inchikeys}
    
    def adapter_status(self) -> dict:
        return {
            "askcos_inventory_adapter": {
                "status": "connected" if self.inventory_df is not None else "not connected",
                "num_compounds": self.count if self.inventory_df is not None else 0,
            }
        }
    
    def inventory_status_details(self, inchikey: str) -> Availability:
        available = self.in_inventory_by_inchikey(inchikey)
        match = self.inventory_df.loc[self.inventory_df[self.inchikey_col] == inchikey]
        source = match[self.source_column].values[0] if self.source_column and not match.empty else None
        url = match[self.url_column].values[0] if self.url_column and not match.empty else None
        return Availability(
            inchikey=inchikey,
            inventory={"available": False},
            commercial_availability={"available": available, "source": source, "url": url}
        )