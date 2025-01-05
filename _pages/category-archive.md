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
anaconda3\envs\[Your Environment Name]\Lib\site-packages\LittleSnowFox
```
Overwrite the files in 
```python
./LittleSnowFox/
```
After installing **LittleSnowFox** using `pip install LittleSnowFox`, the Matlab program will be successfully downloaded. However, the package does not include the accompanying data files. To avoid potential errors in downloading data files through the Python component, you need to manually download the data files. You can directly download the entire `database` folder.

If the folder is placed correctly, you will be prompted with "Do you want to overwrite the `database` folder?" Select "Yes" to write the data files. Ensure that the data directory is correctly overwritten in the Anaconda environment during the process.
  
  
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
  

