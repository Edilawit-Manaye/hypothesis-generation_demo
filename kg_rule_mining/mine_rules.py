# mine_rules_safe.py
from clause import Options, Learner
import os

path_train = "data/kg_triples_light.tsv"
path_rules_out = "output/mined_rules_safe.rules"

os.makedirs(os.path.dirname(path_rules_out), exist_ok=True)

opts = Options()
opts.set("learner.mode", "amie")


opts.set("learner.amie.raw.mins", 5)
opts.set("learner.amie.raw.minhc", 0.0)
opts.set("learner.amie.raw.minpca", 0.0)
opts.set("learner.amie.raw.maxad", 2)
opts.set("learner.amie.raw.const", "")
opts.set("learner.amie.raw.nc", 1)  # use 1 core instead of -threads

learner = Learner(options=opts.get("learner"))

print("Starting AMIE rule mining (safe mode)...")
learner.learn_rules(path_data=path_train, path_output=path_rules_out)

print("Done! Rules written to:", path_rules_out)
