Codes for the deterministic and statistical control of kirigami topology

Reference:
S. Chen, G. P. T. Choi, and L. Mahadevan, 
"Deterministic and statistical control of kirigami topology."
Proceedings of the National Academy of Sciences USA, 2020.

Copyright (c) 2020, Siheng Chen, Gary Pui-Tung Choi, L. Mahadevan

===============================================================
MRP_quad_2x2_all.m: Construct all minimum rigidifying link patterns (MRPs) for a 2x2 quad kirigami
MRP_quad_3x3/4x4/5x5/7x7/3x5.m: Construct a minimum rigidifying link pattern (MRP) for an mxn quad kirigami
MRP_quad_30x30.m: Construct a minimum rigidifying link pattern (MRP) for a 30x30 quad kirigami by hierarchical construction
MRP_quad_3x3_all.m: Construct all minimum rigidifying link patterns (MRPs) for a 3x3 quad kirigami
MRP_quad_3x3_all_with_bdy.m: Construct all minimum rigidifying link patterns (MRPs) for a 3x3 quad kirigami assuming all boundary links are included
MCP_quad_2x2/3x3/4x4.m: Construct a minimum connecting link pattern (MCP) for an mxn quad kirigami
MRP_kagome_5x5/7x7/3x5/5x3.m: Construct a minimum rigidifying link pattern (MRP) for an mxn kagome kirigami
MRP_kagome_all_with_bdy: Construct all minimum rigidifying link patterns (MRPs) for a kagome kirigami assuming all boundary links are included
DoF_NCC_vs_link_density_quad.m: Calculate the degree of freedom (DoF) and the number of connected components (NCC) for a quad kirigami. The number of links added vary from 0 to n_max_link (4*L*(L-1)). The size of the largest connected component and the number of internal rotational DoF can also be calculated.
DoF_NCC_vs_link_density_kagome.m: Similar to calculations above, for the kagome kirigami.
calc_rank.m: Calculate the rank of a sparse matrix using QR method.