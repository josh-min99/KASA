# data_inventory

생성: Stage 0 (`src/s0_download.py`). 총 27건 시도.

| 스터디 | 조직 | 패러다임 | n | 플랫폼 | 희생시각 기록 | 발현행렬 | ISA | 비고 |
|---|---|---|---|---|---|---|---|---|
| OSD-21 | gastrocnemius | ANCHOR | 23 | microarray | X | O | O |  |
| OSD-237 | dorsal skin | HLU | 21 | RNA-seq | X | O | O |  |
| OSD-238 | dorsal skin | FLIGHT | 24 | RNA-seq | X | O | O | 용량초과 건너뜀: GLDS-238_rna_seq_differential_expression_GLbulkRNAseq.csv (134MB) |
| OSD-240 | dorsal skin | FLIGHT | 20 | RNA-seq | X | O | O |  |
| OSD-243 | dorsal skin | FLIGHT | 53 | RNA-seq | X | O | O |  |
| OSD-254 | dorsal skin | FLIGHT | 80 | RNA-seq | X | O | O | 용량초과 건너뜀: GLDS-254_rna_seq_differential_expression_GLbulkRNAseq.csv (446MB) |
| OSD-203 | retina | HLU | 59 | RNA-seq | X | O | O | 용량초과 건너뜀: GLDS-203_rna_seq_differential_expression_GLbulkRNAseq.csv (298MB) |
| OSD-87 | retina | FLIGHT | 6 | microarray | X | O | O |  |
| OSD-194 | retina | FLIGHT | 13 | RNA-seq | X | O | O |  |
| OSD-255 | retina | FLIGHT | 16 | RNA-seq | X | O | O |  |
| OSD-758 | retina | FLIGHT | 35 | RNA-seq | X | O | O |  |
| OSD-876 | gastrocnemius | HLU | 9 | RNA-seq | X | O | O |  |
| OSD-880 | gastrocnemius | HLU | 18 | RNA-seq | X | O | O |  |
| OSD-101 | gastrocnemius | FLIGHT | 16 | RNA-seq | X | O | O |  |
| OSD-419 | gastrocnemius | FLIGHT | 12 | RNA-seq | X | O | O |  |
| OSD-935 | soleus | HLU | 20 | RNA-seq | X | O | O |  |
| OSD-949 | soleus | HLU | 30 | RNA-seq | X | O | O |  |
| OSD-104 | soleus | FLIGHT | 12 | RNA-seq | X | O | O |  |
| OSD-714 | soleus | FLIGHT | 21 | RNA-seq | X | O | O |  |
| OSD-770 | soleus | FLIGHT | 30 | RNA-seq | O | O | O |  |
| OSD-201 | spleen | HLU | 32 | RNA-seq | X | O | O |  |
| OSD-211 | spleen | HLU | 21 | RNA-seq | X | O | O |  |
| OSD-246 | spleen | FLIGHT | 44 | RNA-seq | X | O | O |  |
| OSD-288 | spleen | FLIGHT | 9 | RNA-seq | X | O | O |  |
| OSD-506 | spleen | FLIGHT | 28 | RNA-seq | O | O | O |  |
| OSD-214 | bone marrow | HLU | 24 | RNA-seq | X | O | O |  |
| OSD-690 | bone marrow | FLIGHT | 24 | RNA-seq | X | O | O |  |

## 군 구성

- **OSD-21** (gastrocnemius): Normally Loaded Control | 12 days | calf muscle | Not Applicable=5; Hindlimb Unloaded and Reloaded | 12 days plus reloading for 3.5hr | calf muscle | Not Applicable=5; Hindlimb Unloaded | 12 days | calf muscle | Not Appl
- **OSD-237** (dorsal skin): Normally Loaded Control | non-irradiated=6; Normally Loaded Control | cobalt-57 gamma radiation=6; Hindlimb Unloaded | cobalt-57 gamma radiation=6; Hindlimb Unloaded | non-irradiated=3
- **OSD-238** (dorsal skin): Space Flight | uG | JAXA Chow with FOS=3; Space Flight | uG | JAXA Chow=3; Space Flight | 1G by centrifugation | JAXA Chow with FOS=3; Space Flight | 1G by centrifugation | JAXA Chow=3; Vivarium Control | 1G on Earth | J
- **OSD-240** (dorsal skin): Ground Control=10; Space Flight=10
- **OSD-243** (dorsal skin): Basal Control | 1 | On Earth | Upon euthanasia=10; Ground Control | ~60 | On Earth | Carcass=10; Ground Control | ~30 | On Earth | Upon euthanasia=9; Space Flight | ~30 | On Earth | Upon euthanasia=9; Space Flight | ~60 
- **OSD-254** (dorsal skin): C57BL/6J | Basal Control | 0=10; C3H/HeJ | Basal Control | 0=10; C57BL/6J | Space Flight | ~75=5; C57BL/6J | Ground Control | ~25=5; C57BL/6J | Ground Control | ~75=5; C57BL/6J | Space Flight | ~25=5; C57BL/6J | Vivarium
- **OSD-203** (retina): 7 | non-irradiated | Normally Loaded Control=6; 7 | cobalt-57 gamma radiation | Hindlimb Unloaded=6; 7 | cobalt-57 gamma radiation | Normally Loaded Control=5; 1 | non-irradiated | Normally Loaded Control=5; 4 | cobalt-5
- **OSD-87** (retina): Ground Control=3; Space Flight=3
- **OSD-194** (retina): Space Flight=5; Basal Control=4; Ground Control=4
- **OSD-255** (retina): Ground Control=8; Space Flight=8
- **OSD-758** (retina): Ground Control | 1G on Earth=12; Space Flight | 0.33G by centrifugation=6; Space Flight | 0.66G by centrifugation=6; Space Flight | 1G by centrifugation=6; Space Flight | uG=5
- **OSD-876** (gastrocnemius): Normally Loaded Control=4; Hindlimb Unloaded and Reloaded=3; Hindlimb Unloaded=2
- **OSD-880** (gastrocnemius): Control | Normally Loaded Control=3; Zfp697 mKO | Normally Loaded Control=3; Control | Hindlimb Unloaded=3; Zfp697 mKO | Hindlimb Unloaded=3; Control | Hindlimb Unloaded and Reloaded=3; Zfp697 mKO | Hindlimb Unloaded and
- **OSD-101** (gastrocnemius): Ground Control=7; Space Flight=7; Not Applicable=2
- **OSD-419** (gastrocnemius): Ground Control=8; Space Flight=4
- **OSD-935** (soleus): Normally Loaded Control | sham irradiation | Non-treated=10; Hindlimb Unloaded | x-ray radiation | corticosterone=10
- **OSD-949** (soleus): Control | Female | Normally Loaded Control=5; AMPKalpha mKO | Male | Normally Loaded Control=5; AMPKalpha mKO | Female | Hindlimb Unloaded=5; Control | Male | Hindlimb Unloaded=4; Control | Female | Hindlimb Unloaded=3; 
- **OSD-104** (soleus): Ground Control=6; Space Flight=6
- **OSD-714** (soleus): Ground Control | 1G on Earth | Mouse Habitat Unit 1 (MHU-1)=3; Space Flight | 1G with centrifugation | Mouse Habitat Unit 1 (MHU-1)=3; Space Flight | uG | Mouse Habitat Unit 1 (MHU-1)=3; Ground Control | 1G on Earth | Mo
- **OSD-770** (soleus): Space Flight=10; Ground Control=10; Vivarium Control=10
- **OSD-201** (spleen): Hindlimb Unloaded | not tetanus toxoid injected | not CpG injected=4; Hindlimb Unloaded | tetanus toxoid injected | not CpG injected=4; Hindlimb Loaded Control | tetanus toxoid injected | not CpG injected=4; Hindlimb Unl
- **OSD-211** (spleen): Normally Loaded Control | non-irradiated=6; Normally Loaded Control | cobalt-57 gamma radiation=6; Hindlimb Unloaded | cobalt-57 gamma radiation=6; Hindlimb Unloaded | non-irradiated=3
- **OSD-246** (spleen): Basal Control | 1 | On Earth | Upon euthanasia=10; Ground Control | ~30 | On Earth | Upon euthanasia=9; Space Flight | ~30 | On Earth | Upon euthanasia=8; Ground Control | ~60 | On Earth | Carcass=7; Space Flight | ~60 |
- **OSD-288** (spleen): Space Flight | 1G by centrifugation=3; Ground Control | 1G on Earth=3; Space Flight | uG=3
- **OSD-506** (spleen): Ground Control=10; Space Flight=9; Vivarium Control=9
- **OSD-214** (bone marrow): Hindlimb Unloaded | not tetanus toxoid injected | not CpG injected=3; Hindlimb Unloaded | not tetanus toxoid injected | CpG injected=3; Hindlimb Unloaded | tetanus toxoid injected | not CpG injected=3; Hindlimb Unloaded 
- **OSD-690** (bone marrow): Wild Type | Ground Control=6; Nrf2KO | Ground Control=6; Wild Type | Space Flight=6; Nrf2KO | Space Flight=6