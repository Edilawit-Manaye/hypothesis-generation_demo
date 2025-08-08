# marimo notebook
import marimo

__generated_with = "0.1.0"
app = marimo.App()

# ╔════ Cell 1: Imports
@app.cell
def __():
    import requests
    import json
    import os
    from pronto import Ontology
    import pandas as pd
    
    return requests, json, os, Ontology, pd

# ╔════ Cell 2: Utility Functions
@app.cell
def __():
    def download_file(url, filename):
        if os.path.exists(filename):
            print(f"{filename} already exists, skipping download")
            return
        print(f"Downloading {filename}...")
        response = requests.get(url)
        response.raise_for_status()
        with open(filename, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {filename}")

    def download_json_file(url, filename):
        if os.path.exists(filename):
            print(f"{filename} already exists, loading from local file")
            with open(filename, 'r') as f:
                data = json.load(f)
            return data
        print(f"Downloading {filename}...")
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4)
        print(f"Downloaded {filename}")
        return data

    return download_file, download_json_file



if __name__ == "__main__":
    app.run()