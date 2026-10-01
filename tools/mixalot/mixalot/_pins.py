"""MIXALOT vendor pins — (source id, sha256) for every vendored engine
byte-file. Vendor bytes are never edited in place; a re-vendor updates this
file in the same write. mixalot.verify() re-hashes the vendor set against
PINS fail-closed before any compute and checks the live SOURCES
advisory-loud.
"""

# vendored-relative-path -> (source id, sha256)
PINS = {
    "scrna/telegraph_evaluate.py": (
        "origin:telegraph_evaluate.py", "946da2a61a590d872321b60120f65eca9b44c81d21b5e350fa6daad27a1e367d"),
    "eco/etienne_evaluate.py": (
        "origin:etienne_evaluate.py", "30c48d06e7f692733c652148f119a0d9f79d73b2c1a803dd3e04b7f1c836be8b"),
    "pilot/lsx_direct.py": (
        "origin:lsx_direct.py", "ef1188cee8f2c965f682de971edee2c6a703dd05454978e5b1e7ce9fe7de483c"),
    "pilot/zseries.py": (
        "origin:zseries.py", "5757c734070bb49eef4dde26b5de74fbf7549533a40f9a82d082f62c262a461e"),
    "pilot/closed_form_1var.py": (
        "origin:closed_form_1var.py", "6b800c7a9995055086aafdf7d5bede7c4940c51137217d51ccf000eacdcf8b8c"),
    "pilot/formula_emitter.py": (
        "origin:formula_emitter.py", "0972ecf3bb0f5b3bad68b37f8ae6b9105ffa7289bea33427b9c143001a107baf"),
    "pilot/formula_emitter_dirichlet.py": (
        "origin:formula_emitter_dirichlet.py",
        "530e6f1460a6cd2f33f6b67910e8358cf702b0fe49672d9b455f3eab7055f026"),
    "pilot/lsx55_exact.py": (
        "origin:lsx55_exact.py", "d76c3cbb0832e29a3c901584b5863231c64d825464ac7795b4e91261d50c48fa"),
    "jeff/bigg.py": (
        "origin:bigg.py", "253c4f374cee5fba213cb576cedba2c053c9d0af8c25dd255034abb191394986"),
    "jeff/bigk2.py": (
        "origin:bigk2.py", "d7f422c812a62d3c1030cafaecfb24c3d724fc35a099f22e7beafeffd11fc04d"),
    "jeff/biggen.py": (
        "origin:biggen.py", "74014dbeac952e91654f9dc58926588d4e0f882b2d26793a86c72b54048044f1"),
    "jeff/f_nu.py": (
        "origin:f_nu.py", "1d4b6a23184800dfe43aee31a8b365d7706af1f23c1f514fe824246cf9dd4171"),
    "jeff/w1_collapsed.py": (
        "origin:w1_collapsed.py", "b2bd270826e7575f3d492f3d7e3ed01f86d8d9cbeca45fd0914e116486b83481"),
    "jeff/estimators.py": (
        "origin:estimators.py", "0335271c7ab410a363f2c7a53b2d55cf10399f827acd08524d3cd21c8ad3c318"),
    "tools/annihilator.py": (
        "origin:annihilator.py", "bb0131f94d80367d374984490a2ae0d429005701c71f7120fa9c96c4901a3dd6"),
    "jeff/w1_brute.py": (
        "origin:w1_brute.py", "c8305f90e31da37ae2be58751a31faa8e834a381c2680d2d09e72f1b3bac6e53"),
    "jeff/bigk.py": (
        "origin:bigk.py", "69f0a76981d9d62be7b2a2f31568914a5238490c37d8f6e31ad106bbea0a1e31"),
    "jeff/w4_blind_gf.py": (
        "origin:w4_blind_gf.py", "5e4870aaa89945e8456ddfafcce80928b62d3f66464d82bec88258d02eb270ad"),
    "jeff/w4_dpm_limit.py": (
        "origin:w4_dpm_limit.py", "0d84fcb80f785e3c21d8ebff4ea251f87a5afe338f1c2dd1048f6cd10b1404e3"),
    "jeff/swap_route.py": (
        "origin:swap_route.py", "a65e24dd5d20f4df3f1312de5dcb75c4430163329ebfe837bc179a3dce4c78fd"),
    "census/seg_v1.py": (
        "origin:seg_v1.py", "9101c6e36639c33c8b4b7b52cbb6d61915698eb0c36fdc2c634dd60ebdb2dfe5"),
    "census/seg_blind.py": (
        "origin:seg_blind.py", "276c90766c23eec0a9246a2cc1e51f2ea8c758f0e4801e31649d27f4fdb468e7"),
    "census/mic_power.py": (
        "origin:mic_power.py", "397659beccfbc134cc90784e09ae9eb4fabc92e9fb039505df777bcf8fc0bebf"),
    "census/production_sweep.py": (
        "origin:production_sweep.py", "1596f9d889e66e6c1891b89401b80f12cf1ade725ce27bdceadcc316cfac722a"),
    "blend/frozen_comp_v1.py": (
        "origin:frozen_comp_v1.py", "93384c74011984894274b9b5c723f9f52f7365e2565c5e0ee2ad795bc35b421f"),
    "blend/frozen_comp_blind.py": (
        "origin:frozen_comp_blind.py", "5b12c142508bcc1c9ab933b58c45b5b7e96df94d1363fbeaa1e19f66cfb42e7c"),
    "eco/pilot/ball_engine.py": (
        "origin:ball_engine.py", "54c06b31cb1e010d70525c499d4fb97f348853594fd460cdd6e88706e6c18c6f"),
    "eco/pilot/phase2_certify.py": (
        "origin:phase2_certify.py", "0b33aeb7ec65d5898ed638373b38ae4e04232a9036bee3c6c5d0f83d5e3855cb"),
    "eco/pilot/phase2_multisample_eqI.py": (
        "origin:phase2_multisample_eqI.py", "3aa5917add556c0bfa1d30d7987b7be2072d195801793fa118b1dd3b6e06f1e1"),
    "eco/pilot/multisample.py": (
        "origin:multisample.py", "54765737fa580583f9f0d2a42889be900f3acdc75cac05a2a55d3972ad4ba616"),
    "eco/pilot/multisample_ball.py": (
        "origin:multisample_ball.py", "4e69ae92ca3fcf3dd9917de997847c38a59de4ef33b27c211306f364d3a8d521"),
    "eco/pilot/hier_kron.py": (
        "origin:hier_kron.py", "c5ee8ff77f9e352306525a34936268bb70359eab2bd55b2ac057b59d597983cf"),
    "eco/pilot/hierarchical3.py": (
        "origin:hierarchical3.py", "e8fa14e00a8778bb9bf60ed03442a1c24844213f614bbe3b774b3ba5508cffc8"),
    "eco/pilot/hier_ball.py": (
        "origin:hier_ball.py", "a338967546e8519604d0b5846c876c43ffa0770a34cf8b177cdd183bf3ce8830"),
    "eco/pilot/hier_eval.py": (
        "origin:hier_eval.py", "2234219568b7c5e19e6a6db7f9587930ebbd9728744782af514c0f04e66ad3d3"),
    "eco/pilot/rf_engine.py": (
        "origin:rf_engine.py", "66bf79fb3750796ddbe0661d1ce8b88c073854b6348343dcb0c7b956a74ec090"),
    "eco/pilot/etienne_oracle.py": (
        "origin:etienne_oracle.py", "3761e46c2b2bb9e3a6932d51145f05b5e24627de47114c184cf9ac55fc336dfa"),
    "eco/pilot/gate_engines.py": (
        "origin:gate_engines.py", "b2e0414c4f625b09d302e8c5a33d0e8e1bf4b4c5dadde6be3860715a2e26568c"),
    "eco/pilot/gate_rf.py": (
        "origin:gate_rf.py", "8e478372d3f3fa2c799e58e75f3a712d87a071dced2e92001fb2fd1387f66e97"),
    "eco/pilot/gate_worked_example.py": (
        "origin:gate_worked_example.py", "2066b8b429c39adefbf7a270b6e19e9a2e23f64c2e56b1f708e23017f784e19d"),
    "eco/pilot/verify_hierarchical3.py": (
        "origin:verify_hierarchical3.py", "f8f482b054b2bc348585ce562678139cde3ae3e5133f5774e8e3594ab9640615"),
    "eco/pilot/verify_recurrence.py": (
        "origin:verify_recurrence.py", "10bd1156d2c6c269bf9e2cd8fa47575ca8b587201583a7898ba3bc0577696323"),
}
