from rdkit.Chem import CanonSmiles, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdChemReactions import ReactionFromSmarts


def apply_reaction(substrate: str, reaction: str) -> str:
    """Applies reaction SMARTS to substrate and returns the product smiles

    Args:
        substrate: a SMILES string
        reaction: a reaction SMARTS string

    Returns:
        The product SMILES string
    """
    rd_substrate = MolFromSmiles(substrate)
    rd_reaction = ReactionFromSmarts(reaction)
    products = rd_reaction.RunReactants([rd_substrate])
    products = set([MolToSmiles(product[0]) for product in products])

    return list(products)[0]


def main() -> None:
    """Convert bottromycin precursor to the final product"""
    synthesis = [
        {
            "mite_entry": "MITE0000048",
            "enzyme": "BotP",
            "r_smarts": "[#6:1]-[#16:2]-[#6:3]-[#6:4]-[#6@H:5](-[#7H2:6])-[#6:7](=[O:8])-[#7:9]-[#6H2:10]-[#6:11](=[O:12])-[#7:13]-1-[#6:14]-[#6:15]-[#6:16]-[#6@H:17]-1-[#6:18]=[O:19]>>[#7H2:9]-[#6:10]-[#6:11](=[O:12])-[#7:13]-1-[#6:14]-[#6:15]-[#6:16]-[#6@H:17]-1-[#6:18]=[O:19]",
            "substrate": "CSCC[C@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CC(N)=O)NC(=O)[C@@H]1CCCN1C(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](C)NC(=O)[C@@H](NC(=O)[C@H](CCSC)NC(=O)[C@H](CS)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H]1CCCN1C(=O)CNC(=O)[C@@H](N)CCSC)C(C)C)C(C)C)C(C)C)[C@@H](C)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)NCC(=O)N[C@@H](C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CC(=O)O)C(=O)NCC(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](C)C(=O)N[C@H](C(=O)N[C@@H](CO)C(=O)O)[C@@H](C)O"
        },
        {
            "mite_entry": "MITE0000049",
            "enzyme": "botRMT1",
            "r_smarts": "[#6:1]-[#16:2]-[#6:3]-[#6:4]-[#6@H:5](-[#7:6]-[#6:7](=[O:8])-[#6@H:9](-[#6:10]-[#16:11])-[#7:12]-[#6:13](=[O:14])-[#6@H:15](-[#6:16]-[#6:17](-[#8:19])=[O:18])-[#7:20]-[#6:21](=[O:22])-[#6@H:23](-[#6:24]-[c:25]1[c:26][c:27][c:28][c:29][c:30]1)-[#7:31]-[#6:32](=[O:33])-[#6@@H:34](-[#7:35]-[#6:36](=[O:37])-[#6@@H:38](-[#7:39]-[#6:40](=[O:41])-[#6@@H:42](-[#7:43]-[#6:44](=[O:45])-[#6@@H:46]-1-[#6:47]-[#6:48]-[#6:49]-[#7:50]-1-[#6:51](=[O:52])-[#6H2:53]-[#7H2:54])-[#6:55](-[#6:65])-[#6:66])-[#6:56](-[#6:57])-[#6:58])-[#6:59](-[#6:60])-[#6:61])-[#6:62]=[O:63]>>[#6:1]-[#16:2]-[#6:3]-[#6:4]-[#6@H:5](-[#7:6]-[#6:7](=[O:8])-[#6@H:9](-[#6:10]-[#16:11])-[#7:12]-[#6:13](=[O:14])-[#6@H:15](-[#6:16]-[#6:17](-[#8:19])=[O:18])-[#7:20]-[#6:21](=[O:22])-[#6@@H:23](-[#7:31]-[#6:32](=[O:33])-[#6@@H:34](-[#7:35]-[#6:36](=[O:37])-[#6@@H:38](-[#7:39]-[#6:40](=[O:41])-[#6@@H:42](-[#7:43]-[#6:44](=[O:45])-[#6@@H:46]-1-[#6:47]-[#6:48]-[#6:49]-[#7:50]-1-[#6:51](=[O:52])-[#6H2:53]-[#7H2:54])-[#6:55](-[#6:65])-[#6:66])-[#6:56](-[#6:57])-[#6:58])-[#6:59](-[#6:60])-[#6:61])-[#6@@H:24](-[#6:67])-[c:25]1[c:26][c:27][c:28][c:29][c:30]1)-[#6:62]=[O:63]"
        },
        {
            "mite_entry": "MITE0000046",
            "enzyme": "botRMT2",
            "r_smarts": "[#6:1]-[#16:2]-[#6:3]-[#6:4]/[#6@@:5](-[#6:64]=[#8:65])-[#7:6]-[#6:7](\[#6@@:9](-[#7:12]-[#6:13](/[#6@@:15](-[#7:20]-[#6:21](/[#6@@:23](-[#7:31]-[#6:32](/[#6@:34](-[#6:61](-[#6:63])-[#6:62])-[#7:35]-[#6:36](\[#6@:38](-[#6:58](-[#6:60])-[#6:59])-[#7:39]-[#6:40](/[#6@:42](-[#6:55](-[#6:57])-[#6:56])-[#7:43]-[#6:44](/[#6@:46]1-[#7:50](-[#6:51](-[#6:53]-[#7:54])=[#8:52])-[#6:49]-[#6:48]-[#6:47]-1)=[#8:45])=[#8:41])=[#8:37])=[#8:33])-[#6:24]-[#6:25]1:[#6:30]:[#6:29]:[#6:28]:[#6:27]:[#6:26]:1)=[#8:22])-[#6:16]-[#6:17](-[#8:19])=[#8:18])=[#8:14])-[#6:10]-[#16:11])=[#8:8]>>[#6:1]-[#16:2]-[#6:3]-[#6:4]\[#6@@:5](-[#6:64]=[#8:65])-[#7:6]-[#6:7](=[#8:8])/[#6@:9](-[#6:10]-[#16:11])-[#7:12]-[#6:13](=[#8:14])\[#6@:15](-[#6:16]-[#6:17](=[#8:18])-[#8:19])-[#7:20]-[#6:21](=[#8:22])\[#6@:23](-[#6:24]-[#6:25]1:[#6:30]:[#6:29]:[#6:28]:[#6:27]:[#6:26]:1)-[#7:31]-[#6:32](=[#8:33])/[#6@@:34](-[#7:35]-[#6:36](=[#8:37])/[#6@@:38](-[#7:39]-[#6:40](=[#8:41])\[#6@@:42](-[#7:43]-[#6:44](=[#8:45])/[#6@@:46]1-[#6:47]-[#6:48]-[#6:49]-[#7:50]-1-[#6:51](=[#8:52])-[#6:53]-[#7:54])-[#6:55](-[#6:56])-[#6:57])-[#6:58](-[#6])(-[#6:59])-[#6:60])-[#6:61](-[#6])(-[#6:63])-[#6:62]"
        },
        {
            "mite_entry": "MITE0000047",
            "enzyme": "botRMT3",
            "r_smarts": "[#6:1]-[#16:2]-[#6:3]-[#6:4]/[#6@@:5](-[#6:64]=[#8:65])-[#7:6]-[#6:7](\[#6@@:9](-[#7:12]-[#6:13](/[#6@@:15](-[#7:20]-[#6:21](/[#6@@:23](-[#7:31]-[#6:32](/[#6@:34](-[#6:61](-[#6:63])-[#6:62])-[#7:35]-[#6:36](\[#6@:38](-[#6:58](-[#6:60])-[#6:59])-[#7:39]-[#6:40](/[#6@:42](-[#6:55](-[#6:57])-[#6:56])-[#7:43]-[#6:44](/[#6@:46]1-[#7:50](-[#6:51](-[#6:53]-[#7:54])=[#8:52])-[#6:49]-[#6:48]-[#6:47]-1)=[#8:45])=[#8:41])=[#8:37])=[#8:33])-[#6:24]-[#6:25]1:[#6:30]:[#6:29]:[#6:28]:[#6:27]:[#6:26]:1)=[#8:22])-[#6:16]-[#6:17](-[#8:19])=[#8:18])=[#8:14])-[#6:10]-[#16:11])=[#8:8]>>[#6:1]-[#16:2]-[#6:3]-[#6:4]/[#6@@:5](-[#6:64]=[#8:65])-[#7:6]-[#6:7](\[#6@@:9](-[#7:12]-[#6:13](/[#6@@:15](-[#7:20]-[#6:21](/[#6@@:23](-[#7:31]-[#6:32](/[#6@:34](-[#6:61](-[#6:63])-[#6:62])-[#7:35]-[#6:36](\[#6@:38](-[#6:58](-[#6:59])-[#6:60])-[#7:39]-[#6:40](/[#6@:42](-[#6:55](-[#6:56])-[#6:57])-[#7:43]-[#6:44](/[#6@:46]1-[#7:50](-[#6:51](-[#6:53]-[#7:54])=[#8:52])-[#6:49]-[#6:48]-[#6@@:47]-1/[#6])=[#8:45])=[#8:41])=[#8:37])=[#8:33])-[#6:24]-[#6:25]1:[#6:30]:[#6:29]:[#6:28]:[#6:27]:[#6:26]:1)=[#8:22])-[#6:16]-[#6:17](=[#8:18])-[#8:19])=[#8:14])-[#6:10]-[#16:11])=[#8:8]"
        },
        {
            "mite_entry": "MITE0000052",
            "enzyme": "BotC",
            "r_smarts": "[#6:127]-[#16:126]-[#6:125]-[#6:124]-[#6@H:123](-[#7:128]-[#6:129](=[O:130])-[#6@H:131](-[#6:132]-[#16:133])-[#7:134]-[#6:135](=[O:136])-[#6@H:137](-[#6:138])-[#7:139]-[#6:140](=[O:141])-[#6@H:142](-[#6:143]-[c:144]1[c:145][c:146][c:147][c:148][c:149]1)-[#7:150]-[#6:151](=[O:152])-[#6@@H:153](-[#7:154]-[#6:155](=[O:156])-[#6@@H:157](-[#7:158]-[#6:159](=[O:160])-[#6:161](-[#6:162])-[#7:163]-[#6:164](=[O:165])-[#6@@H:166]-1-[#6:167]-[#6:168]-[#6:169]-[#7:170]-1-[#6:171](=[O:172])-[#6:173]-[#7:174])-[#6:175](-[#6:176])-[#6:177])-[#6:178](-[#6:179])-[#6:180])-[#6:121](=[O:122])-[#7:120]-[#6@@H:119](-[#6@@H:181](-[#6:182])-[#8:183])-[#6:117](=[O:118])-[#7:116]-[#6@@H:114](-[#6:115])-[#6:112](=[O:113])-[#7:111]-[#6@@H:106](-[#6:107]-[#6:108](-[#8:110])=[O:109])-[#6:104](=[O:105])-[#7:103]-[#6@@H:95](-[#6:96]-[c:97]1[c:98][c:99][c:100][c:101][c:102]1)-[#6:93](=[O:94])-[#7:92]-[#6@@H:87](-[#6:88]-[#6:89](-[#6:90])-[#6:91])-[#6:85](=[O:86])-[#7:84]-[#6@@H:82](-[#6:83])-[#6:80](=[O:81])-[#7:79]-[#6@@H:77](-[#6:78])-[#6:75](=[O:76])-[#7:74]-[#6@@H:72](-[#6:73])-[#6:70](=[O:71])-[#7:69]-1-[#6:68]-[#6:67]-[#6:66]-[#6@H:65]-1-[#6:63](=[O:64])-[#7:62]-[#6@@H:57](-[#6:58]-[#6:59](-[#7:60])=[O:61])-[#6:55](=[O:56])-[#7:54]-[#6@@H:49](-[#6:50]-[#6:51](-[#7:52])=[O:53])-[#6:47](=[O:48])-[#7:46]-[#6@@H:44](-[#6:45])-[#6:42](=[O:43])-[#7:41]-[#6@@H:39](-[#6:40])-[#6:37](=[O:38])-[#7:36]-[#6@@H:31](-[#6:32]-[#6:33](-[#6:34])-[#6:35])-[#6:29](=[O:30])-[#7:28]-[#6@@H:25](-[#6:26]-[#8:27])-[#6:23](=[O:24])-[#7:22]-[#6@@H:20](-[#6:21])-[#6:18](=[O:19])-[#7:17]-[#6@@H:15](-[#6:16])-[#6:13](=[O:14])-[#7:12]-[#6@@H:6](-[#6:7]-[#6:8]-[#6:9](-[#8:11])=[O:10])-[#6:4](=[O:5])-[#7:3]-[#6@@H:2](-[#6:1])-[#6:184](=[O:185])-[#7:186]-[#6@@H:187](-[#6:188])-[#6:189](=[O:190])-[#7:191]-[#6@@H:192](-[#6:193])-[#6:194](=[O:195])-[#7:196]-[#6@@H:197](-[#6:198]-[#6:199](-[#6:200])-[#6:201])-[#6:202](=[O:203])-[#7:204]-[#6@@H:205](-[#6:206])-[#6:207](=[O:208])-[#7:209]-[#6@@H:210](-[#6:211]-[#8:212])-[#6:213](=[O:214])-[#7:215]-[#6@@H:216](-[#6:217]-[c:218]1[c:219][n:220][c:221]2[c:222][c:223][c:224][c:225][c:226]12)-[#6:227](=[O:228])-[#7:229]-[#6:230]-[#6:231]=[O:232]>>[#6:127]-[#16:126]-[#6:125]-[#6:124]-[#6@H:123](-[#7:128]-[#6:129](=[O:130])-[#6@@H:131]-1-[#6:132]-[#16:133]-[#6:135](=[#7:134]-1)-[#6@H:137](-[#6:138])-[#7:139]-[#6:140](=[O:141])-[#6@H:142](-[#6:143]-[c:144]1[c:145][c:146][c:147][c:148][c:149]1)-[#7:150]-[#6:151](=[O:152])-[#6@@H:153](-[#7:154]-[#6:155](=[O:156])-[#6@@H:157](-[#7:158]-[#6:159](=[O:160])-[#6:161](-[#6:162])-[#7:163]-[#6:164](=[O:165])-[#6@@H:166]-1-[#6:167]-[#6:168]-[#6:169]-[#7:170]-1-[#6:171](=[O:172])-[#6:173]-[#7:174])-[#6:175](-[#6:176])-[#6:177])-[#6:178](-[#6:179])-[#6:180])-[#6:121](=[O:122])-[#7:120]-[#6@@H:119](-[#6@@H:181](-[#6:182])-[#8:183])-[#6:117](=[O:118])-[#7:116]-[#6@@H:114](-[#6:115])-[#6:112](=[O:113])-[#7:111]-[#6@@H:106](-[#6:107]-[#6:108](-[#8:110])=[O:109])-[#6:104](=[O:105])-[#7:103]-[#6@@H:95](-[#6:96]-[c:97]1[c:98][c:99][c:100][c:101][c:102]1)-[#6:93](=[O:94])-[#7:92]-[#6@@H:87](-[#6:88]-[#6:89](-[#6:90])-[#6:91])-[#6:85](=[O:86])-[#7:84]-[#6@@H:82](-[#6:83])-[#6:80](=[O:81])-[#7:79]-[#6@@H:77](-[#6:78])-[#6:75](=[O:76])-[#7:74]-[#6@@H:72](-[#6:73])-[#6:70](=[O:71])-[#7:69]-1-[#6:68]-[#6:67]-[#6:66]-[#6@H:65]-1-[#6:63](=[O:64])-[#7:62]-[#6@@H:57](-[#6:58]-[#6:59](-[#7:60])=[O:61])-[#6:55](=[O:56])-[#7:54]-[#6@@H:49](-[#6:50]-[#6:51](-[#7:52])=[O:53])-[#6:47](=[O:48])-[#7:46]-[#6@@H:44](-[#6:45])-[#6:42](=[O:43])-[#7:41]-[#6@@H:39](-[#6:40])-[#6:37](=[O:38])-[#7:36]-[#6@@H:31](-[#6:32]-[#6:33](-[#6:34])-[#6:35])-[#6:29](=[O:30])-[#7:28]-[#6@@H:25](-[#6:26]-[#8:27])-[#6:23](=[O:24])-[#7:22]-[#6@@H:20](-[#6:21])-[#6:18](=[O:19])-[#7:17]-[#6@@H:15](-[#6:16])-[#6:13](=[O:14])-[#7:12]-[#6@@H:6](-[#6:7]-[#6:8]-[#6:9](-[#8:11])=[O:10])-[#6:4](=[O:5])-[#7:3]-[#6@@H:2](-[#6:1])-[#6:184](=[O:185])-[#7:186]-[#6@@H:187](-[#6:188])-[#6:189](=[O:190])-[#7:191]-[#6@@H:192](-[#6:193])-[#6:194](=[O:195])-[#7:196]-[#6@@H:197](-[#6:198]-[#6:199](-[#6:200])-[#6:201])-[#6:202](=[O:203])-[#7:204]-[#6@@H:205](-[#6:206])-[#6:207](=[O:208])-[#7:209]-[#6@@H:210](-[#6:211]-[#8:212])-[#6:213](=[O:214])-[#7:215]-[#6@@H:216](-[#6:217]-[c:218]1[c:219][n:220][c:221]2[c:222][c:223][c:224][c:225][c:226]12)-[#6:227](=[O:228])-[#7:229]-[#6:230]-[#6:231]=[O:232]"
        },
        {
            "mite_entry": "MITE0000053",
            "enzyme": "BotCD",
            "r_smarts": "[#6:123]-[#16:122]-[#6:121]-[#6:120]-[#6@H:119](-[#7:124]-[#6:125](=[O:126])-[#6@@H:127]-1-[#6:128]-[#16:131]-[#6:130](=[#7:129]-1)-[#6@H:132](-[#6:133])-[#7:134]-[#6:135](=[O:136])-[#6@H:137](-[#6:138]-[c:139]1[c:140][c:141][c:142][c:143][c:144]1)-[#7:145]-[#6:146](=[O:147])-[#6@H:148](-[#6:171])-[#7:149]-[#6:150](=[O:151])-[#6@H:152](-[#6:170])-[#7:153]-[#6:154](=[O:155])-[#6@H:156](-[#6:169])-[#7:157]-[#6:158](=[O:159])-[#6@@H:160]-1-[#6:161]-[#6:162]-[#6:163]-[#7:164]-1-[#6:165](=[O:166])-[#6:167]-[#7H2:168])-[#6:117](=[O:118])-[#7:116]-[#6@@H:115](-[#6@@H:172](-[#6:173])-[#8:174])-[#6:113](=[O:114])-[#7:112]-[#6@@H:110](-[#6:111])-[#6:108](=[O:109])-[#7:107]-[#6@@H:102](-[#6:103]-[#6:104](-[#8:106])=[O:105])-[#6:100](=[O:101])-[#7:99]-[#6@@H:91](-[#6:92]-[c:93]1[c:94][c:95][c:96][c:97][c:98]1)-[#6:89](=[O:90])-[#7:88]-[#6@@H:86](-[#6:87])-[#6:84](=[O:85])-[#7:83]-[#6@@H:81](-[#6:82])-[#6:79](=[O:80])-[#7:78]-[#6@@H:76](-[#6:77])-[#6:74](=[O:75])-[#7:73]-[#6@@H:68](-[#6:69]-[#6:70](-[#8:72])=[O:71])-[#6:66](=[O:67])-[#7:65]-1-[#6:64]-[#6:63]-[#6:62]-[#6@H:61]-1-[#6:59](=[O:60])-[#7:58]-[#6@@H:53](-[#6:54]-[#6:55](-[#7:56])=[O:57])-[#6:51](=[O:52])-[#7:50]-[#6@@H:45](-[#6:46]-[#6:47](-[#7:48])=[O:49])-[#6:43](=[O:44])-[#7:42]-[#6@@H:40](-[#6:41])-[#6:38](=[O:39])-[#7:37]-[#6@@H:35](-[#6:36])-[#6:33](=[O:34])-[#7:32]-[#6@@H:30](-[#6:31])-[#6:28](=[O:29])-[#7:27]-[#6@@H:25](-[#6:26])-[#6:23](=[O:24])-[#7:22]-[#6@@H:20](-[#6:21])-[#6:18](=[O:19])-[#7:17]-[#6@@H:15](-[#6:16])-[#6:13](=[O:14])-[#7:12]-[#6@@H:6](-[#6:7]-[#6:8]-[#6:9](-[#8:11])=[O:10])-[#6:4](=[O:5])-[#7:3]-[#6@@H:2](-[#6:1])-[#6:175](=[O:176])-[#7:177]-[#6@@H:178](-[#6:179])-[#6:180](=[O:181])-[#7:182]-[#6@@H:183](-[#6:184])-[#6:185](=[O:186])-[#7:187]-[#6@@H:188](-[#6:189]-[#6:190](-[#6:191])-[#6:192])-[#6:193](=[O:194])-[#7:195]-[#6@@H:196](-[#6:197])-[#6:198](=[O:199])-[#7:200]-[#6@@H:201](-[#6:202]-[#8:203])-[#6:204](=[O:205])-[#7:206]-[#6@@H:207](-[#6:208]-[c:209]1[c:210][n:211][c:212]2[c:213][c:214][c:215][c:216][c:217]12)-[#6:218](=[O:219])-[#7:220]-[#6:221]-[#6:222]=[O:223]>>[#6:123]-[#16:122]-[#6:121]-[#6:120]-[#6@H:119](-[#7:124]-[#6:125](=[O:126])-[#6@@H:127]-1-[#6:128]-[#16:131]-[#6:130](=[#7:129]-1)-[#6@H:132](-[#6:133])-[#7:134]-[#6:135](=[O:136])-[#6@H:137](-[#6:138]-[c:139]1[c:140][c:141][c:142][c:143][c:144]1)-[#7:145]-[#6:146](=[O:147])-[#6@H:148](-[#6:171])\[#7:149]=[#6:150]-1/[#7:168]-[#6:167]-[#6:165](=[O:166])-[#7:164]-2-[#6:163]-[#6:162]-[#6:161]-[#6@H:160]-2-[#6:158](=[O:159])-[#7:157]-[#6@@H:156](-[#6:169])-[#6:154](=[O:155])-[#7:153]-[#6@H:152]-1-[#6:170])-[#6:117](=[O:118])-[#7:116]-[#6@@H:115](-[#6@@H:172](-[#6:173])-[#8:174])-[#6:113](=[O:114])-[#7:112]-[#6@@H:110](-[#6:111])-[#6:108](=[O:109])-[#7:107]-[#6@@H:102](-[#6:103]-[#6:104](-[#8:106])=[O:105])-[#6:100](=[O:101])-[#7:99]-[#6@@H:91](-[#6:92]-[c:93]1[c:94][c:95][c:96][c:97][c:98]1)-[#6:89](=[O:90])-[#7:88]-[#6@@H:86](-[#6:87])-[#6:84](=[O:85])-[#7:83]-[#6@@H:81](-[#6:82])-[#6:79](=[O:80])-[#7:78]-[#6@@H:76](-[#6:77])-[#6:74](=[O:75])-[#7:73]-[#6@@H:68](-[#6:69]-[#6:70](-[#8:72])=[O:71])-[#6:66](=[O:67])-[#7:65]-1-[#6:64]-[#6:63]-[#6:62]-[#6@H:61]-1-[#6:59](=[O:60])-[#7:58]-[#6@@H:53](-[#6:54]-[#6:55](-[#7:56])=[O:57])-[#6:51](=[O:52])-[#7:50]-[#6@@H:45](-[#6:46]-[#6:47](-[#7:48])=[O:49])-[#6:43](=[O:44])-[#7:42]-[#6@@H:40](-[#6:41])-[#6:38](=[O:39])-[#7:37]-[#6@@H:35](-[#6:36])-[#6:33](=[O:34])-[#7:32]-[#6@@H:30](-[#6:31])-[#6:28](=[O:29])-[#7:27]-[#6@@H:25](-[#6:26])-[#6:23](=[O:24])-[#7:22]-[#6@@H:20](-[#6:21])-[#6:18](=[O:19])-[#7:17]-[#6@@H:15](-[#6:16])-[#6:13](=[O:14])-[#7:12]-[#6@@H:6](-[#6:7]-[#6:8]-[#6:9](-[#8:11])=[O:10])-[#6:4](=[O:5])-[#7:3]-[#6@@H:2](-[#6:1])-[#6:175](=[O:176])-[#7:177]-[#6@@H:178](-[#6:179])-[#6:180](=[O:181])-[#7:182]-[#6@@H:183](-[#6:184])-[#6:185](=[O:186])-[#7:187]-[#6@@H:188](-[#6:189]-[#6:190](-[#6:191])-[#6:192])-[#6:193](=[O:194])-[#7:195]-[#6@@H:196](-[#6:197])-[#6:198](=[O:199])-[#7:200]-[#6@@H:201](-[#6:202]-[#8:203])-[#6:204](=[O:205])-[#7:206]-[#6@@H:207](-[#6:208]-[c:209]1[c:210][n:211][c:212]2[c:213][c:214][c:215][c:216][c:217]12)-[#6:218](=[O:219])-[#7:220]-[#6:221]-[#6:222]=[O:223]"
        },
        {
            "mite_entry": "MITE0000050",
            "enzyme": "BotAH",
            "r_smarts": "[#6:1]-[#6:2](-[#6:3])-[#6@@H:4]1-[#7:5]-[#6:6](=[O:7])-[#6@H:8](-[#6:9])-[#7:10]-[#6:11](=[O:12])-[#6@H:13](-[#6:14])-[#7:15]-[#6:16](=[O:17])-[#6:18]-[#7:19]\[#6:20]-1=[#7:21]/[#6@@H:22](-[#6:23])-[#6:24](=[O:25])-[#7:26]-[#6@@H:27](-[#6:28])-[#6:29](=[O:30])-[#7:31]-[#6@@H:32](-[#6:33])-[#6:34]-1=[#7:35]-[#6@@H:36](-[#6:37]-[#16:38]-1)-[#6:39](=[O:40])-[#7:41]-[#6@@H:42](-[#6:43])-[#6:44](=[O:45])-[#7:46]-[#6@@H:47](-[#6:48])-[#6:49](=[O:50])-[#7:51]-[#6@@H:52](-[#6:53])-[#6:54](=[O:55])-[#7:56]-[#6@@H:57](-[#6:58])-[#6:59]=[O:60]>>[#6:1]-[#6:2](-[#6:3])-[#6@@H:4]1-[#7:5]-[#6:6](=[O:7])-[#6@H:8](-[#6:9])-[#7:10]-[#6:11](=[O:12])-[#6@H:13](-[#6:14])-[#7:15]-[#6:16](=[O:17])-[#6:18]-[#7:19]\[#6:20]-1=[#7:21]/[#6@@H:22](-[#6:23])-[#6:24](=[O:25])-[#7:26]-[#6@@H:27](-[#6:28])-[#6:29](=[O:30])-[#7:31]-[#6@@H:32](-[#6:33])-[#6:34]-1=[#7:35]-[#6@@H:36](-[#6:37]-[#16:38]-1)-[#6:39](-[#8:40])=[O:41]"
        },
        {
            "mite_entry": "MITE0000051",
            "enzyme": "BotH",
            "r_smarts": "[#6:1]-[#6:2](/[#6@:4]1-[#6:55](=[#8:56])-[#7:54]-[#6@@:53](/[#6:57](-[#6:60])(-[#6:59])-[#6:58])-[#6:18](=[#7:19]\[#6@@:20](-[#6:49](-[#6:52])(-[#6:51])-[#6:50])-[#6:21](-[#7:23]\[#6@:24](-[#6:33](-[#7:35]/[#6@:36](-[#6:41]2-[#16:45]-[#6:44]-[#6@@:43](/[#6:46](=[#8:48])-[#8:47])-[#7:42]=2)-[#6:37]-[#6:38](=[#8:40])-[#8:39])=[#8:34])-[#6@:25](-[#6:27]2:[#6:32]:[#6:31]:[#6:30]:[#6:29]:[#6:28]:2)\[#6:26])=[#8:22])-[#7:17]-[#6:16]-[#6:14](=[#8:15])-[#7:13]2-[#6@@:8](-[#6@@:9](-[#6:11]-[#6:12]-2)/[#6:10])/[#6:6](=[#8:7])-[#7:5]-1)-[#6:3]>>[#6:3]-[#6:2](/[#6@:4]1-[#6:55](=[#8:56])-[#7:54]-[#6@@:53](/[#6:57](-[#6:60])(-[#6:59])-[#6:58])-[#6:18](=[#7:19]\[#6@@:20](-[#6:49](-[#6:51])(-[#6:50])-[#6:52])-[#6:21](-[#7:23]\[#6@:24](-[#6:33](-[#7:35]\[#6@@:36](-[#6:41]2-[#16:45]-[#6:44]-[#6@@:43](/[#6:46](=[#8:48])-[#8:47])-[#7:42]=2)-[#6:37]-[#6:38](=[#8:40])-[#8:39])=[#8:34])-[#6@:25](-[#6:27]2:[#6:28]:[#6:29]:[#6:30]:[#6:31]:[#6:32]:2)\[#6:26])=[#8:22])-[#7:17]-[#6:16]-[#6:14](=[#8:15])-[#7:13]2-[#6@@:8](-[#6@@:9](-[#6:11]-[#6:12]-2)/[#6:10])/[#6:6](=[#8:7])-[#7:5]-1)-[#6:1]"
        },
        {
            "mite_entry": "MITE0000044",
            "enzyme": "BotCYP",
            "r_smarts": "[#6:47]-[#6@H:18](\[#7:17]=[#6:16]-1/[#7:15]-[#6:14]-[#6:12](=[O:13])-[#7:11]-2-[#6:10]-[#6:9]-[#6@@H:7](-[#6:8])-[#6@H:6]-2-[#6:4](=[O:5])-[#7:3]-[#6@@H:2](-[#6:1])-[#6:50](=[O:51])-[#7:49]-[#6@H:48]-1-[#6:52])-[#6:19](=[O:20])-[#7:21]-[#6@@H:22](-[#6@@H:23](-[#6:24])-[c:25]1[c:26][c:27][c:28][c:29][c:30]1)-[#6:31](=[O:32])-[#7:33]-[#6@H:34](-[#6:35]-[#6:36](-[#8:37])=[O:38])-[#6:39]-1=[#7:40]-[#6@@H:41](-[#6:42]-[#16:43]-1)-[#6:44](-[#8H:45])=[O:46]>>[#6:47]-[#6@H:18](\[#7:17]=[#6:16]-1/[#7:15]-[#6:14]-[#6:12](=[O:13])-[#7:11]-2-[#6:10]-[#6:9]-[#6@@H:7](-[#6:8])-[#6@H:6]-2-[#6:4](=[O:5])-[#7:3]-[#6@@H:2](-[#6:1])-[#6:50](=[O:51])-[#7:49]-[#6@H:48]-1-[#6:52])-[#6:19](=[O:20])-[#7:21]-[#6@@H:22](-[#6@@H:23](-[#6:24])-[c:25]1[c:26][c:27][c:28][c:29][c:30]1)-[#6:31](=[O:32])-[#7:33]-[#6@H:34](-[#6:35]-[#6:36](-[#8:37])=[O:38])-[c:39]1[n:40][c:41][c:42][s:43]1"
        },
        {
            "mite_entry": "MITE0000045",
            "enzyme": "BotOMT",
            "r_smarts": "[#6:44]-[#6@H:18](\[#7:17]=[#6:16]-1/[#7:15]-[#6:14]-[#6:12](=[O:13])-[#7:11]-2-[#6:10]-[#6:9]-[#6@@H:7](-[#6:8])-[#6@H:6]-2-[#6:4](=[O:5])-[#7:3]-[#6@@H:2](-[#6:1])-[#6:48](=[O:49])-[#7:47]-[#6@H:45]-1-[#6:46])-[#6:19](=[O:20])-[#7:21]-[#6@@H:22](-[#6@@H:36](-[#6:37])-[c:38]1[c:39][c:40][c:41][c:42][c:43]1)-[#6:23](=[O:24])-[#7:25]-[#6@H:26](-[#6:27]-[#6:28](-[#8H:30])=[O:29])-[c:31]1[n:32][c:33][c:34][s:35]1>>[#6]-[#8:30]-[#6:28](=[O:29])-[#6:27]-[#6@@H:26](-[#7:25]-[#6:23](=[O:24])-[#6@@H:22](-[#7:21]-[#6:19](=[O:20])-[#6@H:18](-[#6:44])\[#7:17]=[#6:16]-1/[#7:15]-[#6:14]-[#6:12](=[O:13])-[#7:11]-2-[#6:10]-[#6:9]-[#6@@H:7](-[#6:8])-[#6@H:6]-2-[#6:4](=[O:5])-[#7:3]-[#6@@H:2](-[#6:1])-[#6:48](=[O:49])-[#7:47]-[#6@H:45]-1-[#6:46])-[#6@@H:36](-[#6:37])-[c:38]1[c:39][c:40][c:41][c:42][c:43]1)-[c:31]1[n:32][c:33][c:34][s:35]1"
        }
    ]
    products = []


    for step in synthesis:
        if len(products) == 0:
            products.append(apply_reaction(substrate=step["substrate"], reaction=step["r_smarts"]))
        else:
            products.append(apply_reaction(substrate=products[-1], reaction=step["r_smarts"]))

    gen_bottromycin_a2 = CanonSmiles(products[-1])
    real_bottromycin_a2 = CanonSmiles("C[C@H]1[C@H]2C(N[C@@H](C(C)C)C(N[C@@H](C(C)(C)C)/C(=N/[C@@H](C(C)(C)C)C(N[C@H](C(N[C@@H](c3sccn3)CC(OC)=O)=O)[C@H](c3ccccc3)C)=O)/NCC(=O)N2CC1)=O)=O")

    print(gen_bottromycin_a2)
    print(real_bottromycin_a2)

    if gen_bottromycin_a2 == real_bottromycin_a2:
        print("True")

if __name__ == "__main__":
    main()
