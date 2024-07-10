import json

from rdkit.Chem import CanonSmiles, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

# ANSI escape codes for colors
RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
RESET = "\033[0m"

class ReactionValidator:
    def __init__(self, json_data):
        self.data = json_data
        self.errors = []
        self.warnings = []

    def validate(self):
        self.validate_json_structure()
        self.validate_reaction_smarts()
        self.validate_validated_reactions()
        return len(self.errors) == 0, self.errors

    def validate_json_structure(self):
        required_keys = ["tailoring", "reactionSMARTS", "validated_reactions"]
        for key in required_keys:
            if key not in self.data.get("reactions", [{}])[0]:
                self.errors.append(f"Missing required key: {key}")

    def validate_reaction_smarts(self):
        reaction_data = self.data.get("reactions", [{}])[0]
        if "reactionSMARTS" not in reaction_data:
            return
        
        reaction_smarts = reaction_data["reactionSMARTS"]["reactionSMARTS"]
        try:
            reaction = ReactionFromSmarts(reaction_smarts)
            if reaction is None:
                self.errors.append("Invalid reaction SMARTS")
        except Exception as e:
            self.errors.append(f"Error parsing reaction SMARTS: {e!s}")

    def validate_validated_reactions(self):
        reaction_data = self.data.get("reactions", [{}])[0]
        if "validated_reactions" not in reaction_data:
            return
        
        for _idx, reaction in enumerate(reaction_data["validated_reactions"]):
            substrate = reaction["substrate_(sub)structure"]
            products = reaction["product_(sub)structure"]
            
            try:
                substrate_mol = MolFromSmiles(self.unescape_smiles(substrate))
                if substrate_mol is None:
                    self.errors.append(f"Invalid substrate {substrate} in validated reaction!")

                expected_mols = set()
                for product in products:
                    product_mol = MolFromSmiles(self.unescape_smiles(product))
                    if product_mol is None:
                        self.errors.append(f"Invalid product {product} in validated reaction!")
                    expected_mols.add(product_mol)

                reaction_smarts = reaction_data["reactionSMARTS"]["reactionSMARTS"]
                reaction = ReactionFromSmarts(reaction_smarts)
                if reaction is None:
                    self.errors.append(f"Invalid reaction SMARTS {reaction_smarts} in validated reaction!")

                predicted_products = set(reaction.RunReactants([substrate_mol]))
                if len(predicted_products) == 0:
                    self.errors.append(f"Reaction did not produce any products in validated reaction!")

                predicted_smiles = set()
                for product in predicted_products:
                    predicted_smiles.add(MolToSmiles(product[0]))

                predicted_mols = set()
                for s in predicted_smiles:
                    predicted_mols.add(MolFromSmiles(self.unescape_smiles(s),sanitize=False))


                expected_smiles = set()
                for mol in expected_mols:
                    expected_smiles.add(self.unescape_smiles(CanonSmiles(MolToSmiles(mol))))

                predicted_smiles = set()
                for mol in predicted_mols:
                    predicted_smiles.add(self.unescape_smiles(CanonSmiles(MolToSmiles(mol))))

                # print(expected_smiles)
                # print(predicted_smiles)

                if not predicted_smiles.intersection(expected_smiles):
                    self.errors.append(f"Predicted products {predicted_smiles} do not match expected products {expected_smiles} in validated reaction!")
            except Exception as e:
                self.warnings.append(f"Exception in validated reaction: {e!s}")

    @staticmethod
    def unescape_smiles(smiles):
        return smiles.replace("\\\\", "\\")


def validate_reaction_json(json_file_path):
    with open(json_file_path) as f:
        json_data = json.load(f)
    
    validator = ReactionValidator(json_data)
    is_valid, errors = validator.validate()
    
    print(f"Validating {json_file_path}:")
    if is_valid:
        print(f"{GREEN}The reaction JSON is valid.{RESET}")
    else:
        print(f"{RED}The reaction JSON is invalid. Errors:{RESET}")
        for error in errors:
            print(f"{RED}- {error}{RESET}")
    
    if validator.warnings:
        print(f"{YELLOW}Warnings:{RESET}")
        for warning in validator.warnings:
            print(f"{YELLOW}- {warning}{RESET}")
    
    print()

# Usage example
if __name__ == "__main__":
    json_files = ["mcjB.json", "mcjC.json", "mibH.json", "rebH.json"]
    for json_file in json_files:
        validate_reaction_json(json_file)
