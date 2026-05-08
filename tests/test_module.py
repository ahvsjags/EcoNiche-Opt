from econiche.module import EcoNicheModule


def test_module_table_has_one_row_per_gene_with_state_and_direction():
    module = EcoNicheModule(
        genes_by_state={
            "tnk_effector": {"GZMB", "PRF1"},
            "caf_ecm_exclusion": {"COL1A1"},
        }
    )
    directions = {"GZMB": -1, "PRF1": -1, "COL1A1": 1}

    table = module.module_table(gene_directions=directions)

    assert set(table["gene"]) == {"GZMB", "PRF1", "COL1A1"}
    assert set(table["state"]) == {"tnk_effector", "caf_ecm_exclusion"}
    assert set(table["direction"]) == {-1, 1}
