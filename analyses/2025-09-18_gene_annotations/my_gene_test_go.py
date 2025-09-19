import mygene
import pandas as pd

# Your full list of gene symbols
gene_symbols = [
    "NOC2L", "WRAP73", "TNFRSF8", "MFAP2", "CAPZB", "DDOST", "ID3", "RPL11", "RSRP1",
    "ATP5IF1", "SNRNP40", "EIF3I", "AKIRIN1", "MACF1", "RPS8", "MUTYH", "UQCRH", "ROR1",
    "SERBP1", "ENSA", "RPS27", "TPM3", "C1orf43", "ADAR", "FLAD1", "ASH1L", "DAP3",
    "ISG20L2", "DHX9-AS1", "NAV1", "KDM5B", "TMEM63A", "PARP1", "PCNX2", "ARID4B",
    "SMYD3", "TP53I3", "PTRHD1", "PPM1G", "SPTBN1", "GFPT1", "MTHFD2", "LRRTM4",
    "TCF7L1", "TGOLN2", "CIAO1", "LIMS1", "TSN", "WDR33", "UBXN4", "ATP5MC3",
    "HSPD1", "BMPR2", "WDR12", "CYP20A1", "METTL21A", "SPAG16", "ITM2C", "NCL",
    "PTMA", "RAB17", "LRRFIP1", "COPS9", "KIF1A", "SEPTIN2", "CRELD1", "RPL32",
    "UBE2E1", "RPL15", "RPL14", "CCDC12", "ELP6", "MAP4", "SHISA5", "PRKAR2A",
    "GPX1", "RBM6", "RPL29", "WDR82", "SPCS1", "TKT", "IL17RD", "PTPRG", "ROBO2",
    "DPPA4", "ITGB5", "RNF7", "EIF4G1", "IGF2BP2", "SENP5", "UBE2K", "RBM47",
    "YTHDC1", "RPL34", "LRBA", "SFRP2", "TMA16", "MSMO1", "TENM3", "ADCY2", "BTF3",
    "HMGCR", "F2R", "TBCA", "EPB41L4A-AS1", "CCDC112", "AP3S1", "SRFBP1", "UQCRQ",
    "FAM13B", "KIF20A", "UBE2D2", "CSNK1A1", "RPS14", "EBF1", "PTTG1", "NPM1",
    "HIGD2A", "PRELID1", "DBN1", "NHP2", "GFPT2", "RACK1", "CDYL", "JARID2", "CAP2",
    "CASC15", "H1-4", "H1-5", "HMGA1", "RPS10", "SNRPC", "RPL10A", "TOMM6",
    "HSP90AB1", "EEF1A1", "PNISR", "HDAC2", "SERINC1", "SF3B5", "UST", "ARID1B",
    "ACTB", "FSCN1", "ZDHHC4", "TRA2A", "FKBP14", "COA1", "PPIA", "CCT6A", "GALNT17",
    "BUD23", "ABHD11", "CLDN3", "EIF4H", "GTF2I", "POR", "STYXL1", "HSPB1",
    "ARPC1B", "PDAP1", "ATP5MF", "ZKSCAN1", "PLOD3", "CUX1", "SRPK2", "ARF5", "CALD1",
    "FMC1", "XRCC2", "UBE3C", "AGPAT5", "XKR6", "PLPP5", "FGFR1", "SFRP1", "PRKDC",
    "CHD7", "FABP5", "ESRP1", "PABPC1", "MTBP", "GFUS", "TONSL", "RPL8", "RPS6",
    "NDUFB6", "PIGO", "CLTA", "TOMM5", "PRPF4", "ATP6V1G1", "NR6A1", "RPL35", "SCAI",
    "RPL12", "BBLN", "SET", "PTPA", "C9orf78", "PRRC2B", "EHMT1", "NUDT5", "YME1L1",
    "PARD3", "JMJD1C", "PSAP", "CAMK2G", "KIF11", "PPRC1", "ACTR1A", "TCF7L2", "RNH1",
    "TALDO1", "RPLP2", "SBF2", "MDK", "MTCH2", "TMX2", "FADS1", "FTH1", "EEF1G",
    "VPS51", "FAU", "DRAP1", "BANF1", "SF3B2", "GSTP1", "RPS3", "PAK1", "CLNS1A",
    "FDX1", "THY1", "GPC6", "COL4A2", "CUL4A", "MRPL52", "PSMB5", "MED6", "FOXN3",
    "MTA1", "SNRPN", "SNHG14", "UBE3A", "HERC2", "ZNF106", "PDIA3", "SERF2", "DUT",
    "EID1", "TPM1", "RPL4", "RPLP1", "MAN2C1", "RPS17", "MFGE8", "AP3S2", "MRPL28",
    "NME4", "METTL26", "RHOT2", "HAGHL", "MRPS34", "FAHD1", "RPS2", "CCNF", "ELOB",
    "ZSCAN10", "CARHSP1", "PAGR1", "NFAT5", "DDX19A", "GLG1", "USP10", "RPL13",
    "PITPNA", "PAFAH1B1", "SLC25A11", "PFN1", "CTDNEP1", "EIF5A", "PLSCR3", "SNHG29",
    "DRG2", "RPL23A", "TAF15", "BRCA1", "SLC25A39", "KIF18B", "CBX1", "IGF2BP1",
    "NME2", "TRIM37", "PITPNC1", "BPTF", "RPL38", "MRPL58", "JPT1", "SAP30BP",
    "CHMP6", "MYL12B", "DLGAP1", "RPL17", "TXNL1", "CNDP2", "BSG", "WDR18",
    "POLR2E", "FAM174C", "UQCR11", "TCF3", "BTBD2", "LSM7", "TIMM13", "EEF2", "RPL36",
    "HNRNPM", "EIF3G", "MRPL4", "ILF3", "PRKCSH", "ZNF791", "GET3", "GADD45GIP1",
    "NDUFB7", "BRD4", "RPL18A", "PIK3R2", "LSM4", "FKBP8", "COPE", "DDX49", "NDUFA13",
    "ZNF714", "ZNF730", "ZNF91", "POLR2I", "TBCB", "LINC00665", "ACTN4", "PLD3",
    "SHKBP1", "RPS19", "APOE", "SNRPD2", "AP2S1", "NAPA", "KDELR1", "RPL18",
    "PIH1D1", "RPL13A", "RPS11", "PRR12", "PRMT1", "PNKP", "PPP2R1A", "ZNF160",
    "MYADM", "CACNG8", "NDUFA3", "PRPF31", "RPS9", "U2AF2", "EPN1", "TRIM28",
    "SNRPB", "MRPS26", "ADISSP", "HM13", "EIF2S2", "RBM39", "MYL9", "YWHAB", "SDC4",
    "PFDN4", "PPDPF", "ITSN1", "HMGN1", "U2AF1", "UBE2G2", "CECR2", "UBE2L3",
    "NIPSNAP1", "PRR14L", "YWHAH", "LARGE1", "TXN2", "ANKRD54", "RPL3", "ST13",
    "RRP7BP", "LINC00685", "TIMM17B", "GNL3L", "PBDC1", "FIRRE", "LDOC1", "RPL10",
    "MT-CO1", "MT-CO2", "MT-CO3"
]

# Initialize MyGeneInfo client
mg = mygene.MyGeneInfo()

# Define prefixes for ribosomal proteins
ribosomal_prefixes = ['RPL', 'RPS', 'MRPL', 'MRPS']

# Query the database
gene_annotations = mg.querymany(
    gene_symbols,
    scopes='symbol',
    fields='go,name',
    species='human'
)

# Process the results
ribo_genes = []
for result in gene_annotations:
    gene_symbol = result.get('query')

    # Heuristic Check: Does the symbol start with a ribosomal prefix?
    is_ribo_by_symbol = any(gene_symbol.startswith(prefix) for prefix in ribosomal_prefixes)

    # Annotation Check: Does the name or GO term contain "ribosomal"?
    is_ribo_by_annotation = False
    
    # Check the full name
    name = result.get('name', '').lower()
    if 'ribosomal' in name:
        is_ribo_by_annotation = True

    # Check the GO terms
    go_data = result.get('go')
    if go_data:
        # The 'go' field can be a dictionary or a list of dictionaries
        if isinstance(go_data, list):
            for go_entry in go_data:
                if 'ribosome' in go_entry.get('term', '').lower():
                    is_ribo_by_annotation = True
                    break
        elif isinstance(go_data, dict):
            if 'ribosome' in go_data.get('term', '').lower():
                is_ribo_by_annotation = True

    # If either check is true, add the gene to the list
    if is_ribo_by_symbol or is_ribo_by_annotation:
        ribo_genes.append(gene_symbol)

print("The following genes are annotated as ribosomal:")
print(ribo_genes)
print(f"Total ribosomal genes found: {len(ribo_genes)}")