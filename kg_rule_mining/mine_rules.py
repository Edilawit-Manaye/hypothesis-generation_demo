from clause import Options, Learner
import os
import sys
import logging
from datetime import datetime

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('output/mining_log.txt')
        ]
    )
    return logging.getLogger(__name__)

def validate_input_file(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input data file not found: {file_path}")
    if os.path.getsize(file_path) == 0:
        raise ValueError(f"Input data file is empty: {file_path}")
    return True

def create_output_directory(dir_path):
    try:
        os.makedirs(dir_path, exist_ok=True)
        return True
    except OSError as e:
        raise OSError(f"Failed to create output directory {dir_path}: {e}")

def configure_amie_options():
    opts = Options()
    opts.set("learner.mode", "amie")
    opts.set("learner.amie.raw.mins", 5)
    opts.set("learner.amie.raw.minhc", 0.0)
    opts.set("learner.amie.raw.minpca", 0.0)
    opts.set("learner.amie.raw.maxad", 2)
    opts.set("learner.amie.raw.const", "")
    opts.set("learner.amie.raw.nc", 1)
    return opts

def initialize_learner(options):
    return Learner(options=options.get("learner"))

def log_mining_parameters(logger, path_train, path_rules_out):
    logger.info("AMIE Rule Mining Parameters:")
    logger.info(f"Input data: {path_train}")
    logger.info(f"Output rules: {path_rules_out}")
    logger.info("Configuration: mins=5, minhc=0.0, minpca=0.0, maxad=2, nc=1")

def main():
    logger = setup_logging()
    
    path_train = "data/kg_triples_light.tsv"
    path_rules_out = "output/mined_rules_safe.rules"
    
    try:
        logger.info("Starting AMIE rule mining process...")
        
        validate_input_file(path_train)
        create_output_directory(os.path.dirname(path_rules_out))
        
        opts = configure_amie_options()
        learner = initialize_learner(opts)
        
        log_mining_parameters(logger, path_train, path_rules_out)
        
        logger.info("Beginning rule mining operation...")
        start_time = datetime.now()
        
        learner.learn_rules(path_data=path_train, path_output=path_rules_out)
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        logger.info(f"Rule mining completed successfully in {duration:.2f} seconds")
        logger.info(f"Output written to: {path_rules_out}")
        
        if os.path.exists(path_rules_out):
            file_size = os.path.getsize(path_rules_out)
            logger.info(f"Output file size: {file_size} bytes")
        
        print("Done! Rules written to:", path_rules_out)
        
    except FileNotFoundError as e:
        logger.error(f"File error: {e}")
        sys.exit(1)
    except OSError as e:
        logger.error(f"OS error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error during rule mining: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
