##### Wildcard constraints #####
#wildcard_constraints:
#    vartype="snvs|indels",
#    sample="|".join(samples.index),
#    unit="|".join(units["unit"]),
#    contig="|".join(contigs)


##### Helper functions #####

def get_resource(rule,resource):
    try:
        return config["resources"][rule][resource]
    except KeyError:
        return config["resources"]["default"][resource]
