---
title: "Install & Guide"
permalink: /categories/
author_profile: true
---



Dependencies
---
Python = `3.9.18`(No compatibility testing performed)  
Matlab = `R2024a`(In theory, all versions are supported, because Matlab is backward compatible)   
R = `4.3.2`  
Anaconda/Miniconda (Highly recommend)  
  
  
PIP Install
---
```python
pip install LittleSnowFox
```  


Download the database
---
```python
https://biocomputing.cowtransfer.com/s/b7f5aa9cc9ee4e  
```  
Password  
```python
fl1n09
```   
When installing TOGGLE (KALIN) using Anaconda, the path is typically:  
```python
./anaconda3/envs/[Your Environment Name]/Lib/site-packages/LittleSnowFox
```
Overwrite the files in 
```python
./LittleSnowFox/
```
After installing **LittleSnowFox** using `pip install LittleSnowFox`, the Matlab program will be successfully downloaded. However, the package does not include the accompanying data files. To avoid potential errors in downloading data files through the Python component, you need to manually download the data files. You can directly download the entire `database` folder.

If the folder is placed correctly, you will be prompted with "Do you want to overwrite the `database` folder?" Select "Yes" to write the data files. Ensure that the data directory is correctly overwritten in the Anaconda environment during the process.

```python
#Database folder Structure：
---[LittleSnowFox]
     ---[database]
     ---[Clustering_sample]
             ---[e-cigarettes], [e-cigarettes_difference_expression], [fibroblasts]
                    ---[data] #File processed by Jupyter Notebook
                    ---[Rdata] #File processed by R
                    ---[result] #File processed by Matlab
                    ---main_v3_[Sample Name].m #Matlab Main Coding
     ---[Tracing_sample] 
            ---[Hematopoiesis], [Nerveferroptosis], [Reprogramming]
                    ---[data] #File processed by Jupyter Notebook
                    ---[Rdata] #File processed by R
                    ---[result] #File processed by Matlab
```   

Cover Letter
---  
Summary:  
This study developed TOGGLE to address challenges in cell fate analysis. Using hematopoietic and fibroblast reprogramming datasets, TOGGLE reconstructed cell lineages, identified progenitor cells, and predicted developmental trajectories without extensive preprocessing. It classified functionally continuous cells into distinct types and tracked intra-class differentiation. Additionally, TOGGLE analyzed RNA pathway changes in e-cigarette and nicotine exposure datasets, highlighting non-genetic and metabolic influences. By identifying RNA pathway transformations and cellular function execution, TOGGLE offers a robust framework for understanding programmed cell fate transitions.  
Meanwhile, we used neurons undergoing ferroptosis as sample to establish the programmed fate lineage of mature cells. We mapped the trajectory of neuronal programmed cell death and identified neuron populations primarily undergoing ferroptosis. Our findings establish a novel framework for investigating ferroptosis(programmed cell fate) mechanisms and offer TOGGLE as a powerful tool for mapping cell death trajectories across tissues. This study lays the groundwork for future therapeutic strategies targeting ferroptosis in ischemic stroke and other diseases involving programmed cell death.  
Compared with previous algorithm, TOGGLE addresses tracing mature cell fates, overcoming challenges in lineage tracing where mature cells lack clear fate boundaries and often exist in highly similar temporal states. It reconstructs complete fate trajectories, grouping cells by fate progression without pre-identifying progenitor or descendant cells. Based on a corollary of Takens' theorem, TOGGLE enables dynamic reconstruction of cell fate, including differentiation lineages (e.g., hematopoietic cells), reprogramming (e.g., induced cells), programmed fates (e.g., ferroptosis), and functional transformations (e.g., e-cigarette-induced cardiac changes). This tool establishes comprehensive fate-tracking for mature cells and their functions.   

Biology Challenges:  
Terminally differentiated cells undergo programmed cell death through distinct molecular pathways, including ferroptosis, apoptosis, and pyroptosis. However, identifying specific types remains challenging due to the complexity of signaling pathways, marker genes, and the limitations of current tools.   
Algorithm Challenges:  
•	Endpoint loss (e.g., missing terminal samples like dead neurons in stroke).  
•	Time loss (difficulty in sampling at multiple intervals).  
•	Sample imbalance (uneven representation of cell stages), and the absence of clear fate boundaries (gradual transitions in processes like ferroptosis).   
•	Lack of Boundaries: cell soft and same type and subtype lack type boundaries, fate prediction of neurons lacks clear boundaries between cell states  
Innovation of Methods:  
1.	Based on the ability to perform mature cell subtype/type classification and lineage tracing prediction, we extended functionality to enable classification within the same subtype/type of cells based on their specific functional states.  
2.	While existing techniques can only conduct subtype classification and lineage tracing, we enhanced the resolution of these principles to enable programmed fate differentiation, distinguishing different stages of cell fate.  
3.	We demonstrated the potential to detect RNA pathway alterations in the e-cigarette dataset and the ability to classify functional cell interactions in the rat myocardial infarction (MI) model dataset.  

Terminology Explanation:  
•	Subtype Classification: A method of categorizing cell types. Specific mature cells can differentiate into different subtypes under various conditions. For instance, CD4⁺ T cells can transition into Th1 cells, Th2 cells, Treg cells, or Tfh cells. These all belong to the broader category of T cells but represent distinct subcategories. For clarity, we define that if a cell type does not generate other subtypes, it is considered to have only one subtype.  
•	Lineage Tracing: A functional classification method for cells. It predicts the types of cells that a progenitor cell with differentiation potential will develop into. For example, hematopoietic cells may differentiate into monocytes or neutrophils.  
•	Programmed Fate Differentiation: A classification of cell fate stages. Within the same cell subtype or type, specific programmed events are initiated for a particular purpose, such as healthy cells undergoing apoptosis or ferroptosis. This is a new topic we trying to talk for single cell.
Discover of Biology:  
•	We mapped the trajectory of neuronal programmed cell death and identified neuron populations primarily undergoing ferroptosis.  
•	Key ferroptosis-driving genes, including Ctsb, Mtdh, Ndrg1, Smad7, and others, were highlighted as potential therapeutic targets. By characterizing these pathways, we provide new insights into the regulation of ferroptosis, which is essential for rescuing neurons and mitigating damage in affected brain areas.  

Advanced Comparison: 

<br>
<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Advanced.png" 
     alt="Advanced.png" 
     title="Advanced.png">



  
Raw data
---
We provide a format-converted dataset on the cloud drive.  
If you wish to access the original dataset, you can do so through the following link.  
𝙏𝙝𝙖𝙣𝙠𝙨 𝙩𝙤 𝙩𝙝𝙚𝙨𝙚 𝙖𝙪𝙩𝙝𝙤𝙧𝙨 𝙛𝙤𝙧 𝙩𝙝𝙚𝙞𝙧 𝙨𝙥𝙞𝙧𝙞𝙩 𝙤𝙛 𝙨𝙝𝙖𝙧𝙞𝙣𝙜. 𝙏𝙝𝙚𝙮 𝙖𝙧𝙚 𝙜𝙧𝙚𝙖𝙩 𝙥𝙚𝙤𝙥𝙡𝙚.    
  
**Hematopoietic Development (GSE140802)**  
𝘽𝙪𝙣𝙣𝙚, 𝘾., 𝙎𝙩𝙖𝙧𝙠, 𝙎.𝙂., 𝙂𝙪𝙩, 𝙂. 𝙚𝙩 𝙖𝙡. 𝙇𝙚𝙖𝙧𝙣𝙞𝙣𝙜 𝙨𝙞𝙣𝙜𝙡𝙚-𝙘𝙚𝙡𝙡 𝙥𝙚𝙧𝙩𝙪𝙧𝙗𝙖𝙩𝙞𝙤𝙣 𝙧𝙚𝙨𝙥𝙤𝙣𝙨𝙚𝙨 𝙪𝙨𝙞𝙣𝙜 𝙣𝙚𝙪𝙧𝙖𝙡 𝙤𝙥𝙩𝙞𝙢𝙖𝙡 𝙩𝙧𝙖𝙣𝙨𝙥𝙤𝙧𝙩. 𝙉𝙖𝙩 𝙈𝙚𝙩𝙝𝙤𝙙𝙨 𝟮𝟬, 𝟭𝟳𝟱𝟵–𝟭𝟳𝟲𝟴 (𝟮𝟬𝟮𝟯).    
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140802  

  
**Fibroblast Reprogramming (GSE99915)**  
𝙉𝙚𝙘𝙯𝙮𝙥𝙤𝙧, 𝙀𝙫𝙖𝙣 𝙒., 𝙚𝙩 𝙖𝙡. "𝙚-𝘾𝙞𝙜𝙖𝙧𝙚𝙩𝙩𝙚 𝙖𝙚𝙧𝙤𝙨𝙤𝙡 𝙧𝙚𝙙𝙪𝙘𝙚𝙨 𝙡𝙚𝙛𝙩 𝙫𝙚𝙣𝙩𝙧𝙞𝙘𝙪𝙡𝙖𝙧 𝙛𝙪𝙣𝙘𝙩𝙞𝙤𝙣 𝙞𝙣 𝙖𝙙𝙤𝙡𝙚𝙨𝙘𝙚𝙣𝙩 𝙢𝙞𝙘𝙚." 𝘾𝙞𝙧𝙘𝙪𝙡𝙖𝙩𝙞𝙤𝙣 𝟭𝟰𝟱.𝟭𝟭 (𝟮𝟬𝟮𝟮): 𝟴𝟲𝟴-𝟴𝟳𝟬.  
𝘽𝙞𝙙𝙙𝙮 𝘽𝘼, 𝙆𝙤𝙣𝙜 𝙒, 𝙆𝙖𝙢𝙞𝙢𝙤𝙩𝙤 𝙆, 𝙂𝙪𝙤 𝘾 𝙚𝙩 𝙖𝙡. 𝙎𝙞𝙣𝙜𝙡𝙚-𝙘𝙚𝙡𝙡 𝙢𝙖𝙥𝙥𝙞𝙣𝙜 𝙤𝙛 𝙡𝙞𝙣𝙚𝙖𝙜𝙚 𝙖𝙣𝙙 𝙞𝙙𝙚𝙣𝙩𝙞𝙩𝙮 𝙞𝙣 𝙙𝙞𝙧𝙚𝙘𝙩 𝙧𝙚𝙥𝙧𝙤𝙜𝙧𝙖𝙢𝙢𝙞𝙣𝙜. 𝙉𝙖𝙩𝙪𝙧𝙚 𝟮𝟬𝟭𝟴 𝘿𝙚𝙘;𝟱𝟲𝟰(𝟳𝟳𝟯𝟱):𝟮𝟭𝟵-𝟮𝟮𝟰. 𝙋𝙈𝙄𝘿: 𝟯𝟬𝟱𝟭𝟴𝟴𝟱𝟳  
𝙆𝙖𝙢𝙞𝙢𝙤𝙩𝙤 𝙆, 𝘼𝙙𝙞𝙡 𝙈𝙏, 𝙅𝙞𝙣𝙙𝙖𝙡 𝙆, 𝙃𝙤𝙛𝙛𝙢𝙖𝙣𝙣 𝘾𝙈 𝙚𝙩 𝙖𝙡. 𝙂𝙚𝙣𝙚 𝙧𝙚𝙜𝙪𝙡𝙖𝙩𝙤𝙧𝙮 𝙣𝙚𝙩𝙬𝙤𝙧𝙠 𝙧𝙚𝙘𝙤𝙣𝙛𝙞𝙜𝙪𝙧𝙖𝙩𝙞𝙤𝙣 𝙞𝙣 𝙙𝙞𝙧𝙚𝙘𝙩 𝙡𝙞𝙣𝙚𝙖𝙜𝙚 𝙧𝙚𝙥𝙧𝙤𝙜𝙧𝙖𝙢𝙢𝙞𝙣𝙜. 𝙎𝙩𝙚𝙢 𝘾𝙚𝙡𝙡 𝙍𝙚𝙥𝙤𝙧𝙩𝙨 𝟮𝟬𝟮𝟯 𝙅𝙖𝙣 𝟭𝟬;𝟭𝟴(𝟭):𝟵𝟳-𝟭𝟭𝟮. 𝙋𝙈𝙄𝘿: 𝟯𝟲𝟱𝟴𝟰𝟲𝟴𝟱     
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99915  

    
**E-Cigarette and Nicotine Exposures (GSE183614)**   
𝙉𝙚𝙘𝙯𝙮𝙥𝙤𝙧 𝙀𝙒, 𝙎𝙖𝙡𝙙𝙖𝙣̃𝙖 𝙏𝘼, 𝙈𝙚𝙖𝙧𝙨 𝙈𝙅, 𝘼𝙨𝙡𝙖𝙣𝙚𝙧 𝘿𝙈 𝙚𝙩 𝙖𝙡. 𝙚-𝘾𝙞𝙜𝙖𝙧𝙚𝙩𝙩𝙚 𝘼𝙚𝙧𝙤𝙨𝙤𝙡 𝙍𝙚𝙙𝙪𝙘𝙚𝙨 𝙇𝙚𝙛𝙩 𝙑𝙚𝙣𝙩𝙧𝙞𝙘𝙪𝙡𝙖𝙧 𝙁𝙪𝙣𝙘𝙩𝙞𝙤𝙣 𝙞𝙣 𝘼𝙙𝙤𝙡𝙚𝙨𝙘𝙚𝙣𝙩 𝙈𝙞𝙘𝙚. 𝘾𝙞𝙧𝙘𝙪𝙡𝙖𝙩𝙞𝙤𝙣 𝟮𝟬𝟮𝟮 𝙈𝙖𝙧 𝟭𝟱;𝟭𝟰𝟱(𝟭𝟭):𝟴𝟲𝟴-𝟴𝟳𝟬. 𝙋𝙈𝙄𝘿: 𝟯𝟱𝟭𝟴𝟰𝟱𝟳𝟬  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183614  
  
    
**Rat Myocardial Infarction (MI) model**  
𝙏𝙪𝙛𝙖𝙣 𝙏, 𝘾𝙤𝙢𝙚𝙧𝙩𝙥𝙖𝙮 𝙂, 𝙑𝙞𝙡𝙡𝙖𝙣𝙞 𝘼, 𝙉𝙚𝙡𝙨𝙤𝙣 𝙂𝙈 𝙚𝙩 𝙖𝙡. 𝙍𝙖𝙥𝙞𝙙 𝙪𝙣𝙡𝙚𝙖𝙨𝙝𝙞𝙣𝙜 𝙤𝙛 𝙢𝙖𝙘𝙧𝙤𝙥𝙝𝙖𝙜𝙚 𝙚𝙛𝙛𝙚𝙧𝙤𝙘𝙮𝙩𝙞𝙘 𝙘𝙖𝙥𝙖𝙘𝙞𝙩𝙮 𝙫𝙞𝙖 𝙩𝙧𝙖𝙣𝙨𝙘𝙧𝙞𝙥𝙩𝙞𝙤𝙣𝙖𝙡 𝙥𝙖𝙪𝙨𝙚 𝙧𝙚𝙡𝙚𝙖𝙨𝙚. 𝙉𝙖𝙩𝙪𝙧𝙚 𝟮𝟬𝟮𝟰 𝘼𝙥𝙧;𝟲𝟮𝟴(𝟴𝟬𝟬𝟳):𝟰𝟬𝟴-𝟰𝟭𝟱. 𝙋𝙈𝙄𝘿: 𝟯𝟴𝟰𝟴𝟬𝟴𝟴𝟯  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253578  

**Programmed Cell Death Trajectories Prediction in Neurons**   
𝙉𝙖𝙠𝙖𝙢𝙪𝙧𝙖, 𝘼𝙠𝙖𝙧𝙞, 𝙚𝙩 𝙖𝙡. "𝙋𝙇𝘼𝟮𝙂𝟮𝙀-𝙢𝙚𝙙𝙞𝙖𝙩𝙚𝙙 𝙡𝙞𝙥𝙞𝙙 𝙢𝙚𝙩𝙖𝙗𝙤𝙡𝙞𝙨𝙢 𝙩𝙧𝙞𝙜𝙜𝙚𝙧𝙨 𝙗𝙧𝙖𝙞𝙣-𝙖𝙪𝙩𝙤𝙣𝙤𝙢𝙤𝙪𝙨 𝙣𝙚𝙪𝙧𝙖𝙡 𝙧𝙚𝙥𝙖𝙞𝙧 𝙖𝙛𝙩𝙚𝙧 𝙞𝙨𝙘𝙝𝙚𝙢𝙞𝙘 𝙨𝙩𝙧𝙤𝙠𝙚." 𝙉𝙚𝙪𝙧𝙤𝙣 𝟭𝟭𝟭.𝟭𝟵 (𝟮𝟬𝟮𝟯): 𝟮𝟵𝟵𝟱-𝟯𝟬𝟭𝟬.  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232429  
  

