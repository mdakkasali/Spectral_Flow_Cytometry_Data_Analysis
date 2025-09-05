# Spectral Flow Cytometry Data Analysis

## Publication Summary

This repository accompanies the paper:
**"Development of a Spectral Flow Cytometry Analysis Pipeline for High-dimensional Immune Cell Characterization"**

Published in *The Journal of Immunology*, December 1, 2024.

## Citation

**Authors:** Donald Vardaman, III; Md Akkas Ali; Md Hasanul Banna Siam; Chase Bolding; Harrison Tidwell; Holly R Stephens; Mallikarjun Patil; Daniel J Tyrrell

**Journal:** The Journal of Immunology

**Year:** 2024

**Date:** December 1, 2024

**DOI:** 10.4049/jimmunol.2400370

**Open Access PDF:** [Link to PDF](https://academic.oup.com/jimmunol/article-pdf/213/11/1713/61496676/ji2400370.pdf)

---

## Code Architecture

### Repository Structure

- **Jupyter Notebooks**: Interactive analysis workflows with step-by-step documentation
  - `4_16_24_spleen_channel-Copy1.ipynb`: Main analysis notebook
  - Additional supporting notebooks for specialized analyses
- **Data Directory**: Organized storage for spectral flow cytometry datasets
- **Utilities**: Helper scripts and functions for data processing
- **requirements.txt**: Complete list of Python dependencies
- **Documentation**: Comprehensive guides and method descriptions

### Key Libraries and Packages

- **Python**: Core programming language (version 3.x)
- **pandas**: Data manipulation and analysis
- **numpy**: Numerical computing and array operations
- **scanpy**: Single-cell analysis toolkit adapted for flow cytometry
- **matplotlib**: Data visualization and plotting
- **seaborn**: Statistical data visualization
- **scikit-learn**: Machine learning algorithms for clustering and dimensionality reduction
- **flowio**: Flow cytometry file reading and processing
- **anndata**: Annotated data matrices for high-dimensional analysis

### Analysis Workflow

1. **Data Loading**: Import spectral flow cytometry files (.fcs format)
2. **Preprocessing**: Quality control, filtering, and initial data cleaning
3. **Batch Correction**: Harmonization across experimental batches
4. **Clustering**: Unsupervised identification of immune cell populations
5. **Statistical Analysis**: Differential expression and population comparisons
6. **Visualization**: Generation of publication-ready figures and plots

---

## Methods

The analysis pipeline reflects the methodological approach detailed in the publication:

### 1. Sample Preparation

- Tissue processing and cell isolation protocols
- Antibody panel design and optimization
- Sample staining and fixation procedures

### 2. Spectral Flow Cytometry

- Multi-parameter spectral flow cytometry data acquisition
- High-dimensional immune cell characterization
- Quality control measures for spectral data

### 3. Data Processing

- Raw data import and initial quality assessment
- Compensation and spectral unmixing
- Event filtering and debris removal

### 4. Quality Control

- Batch effect identification and assessment
- Sample quality metrics evaluation
- Data integrity validation

### 5. Clustering Analysis

- Unsupervised clustering algorithms (e.g., Leiden, Louvain)
- Dimensionality reduction techniques (UMAP, t-SNE)
- Cell population identification and annotation

### 6. Statistical Analysis

- Differential abundance testing
- Marker expression analysis
- Statistical significance assessment
- Multiple comparison corrections

---

## Usage Instructions

### Prerequisites

- Python 3.8 or higher
- Jupyter Notebook or JupyterLab environment
- Sufficient RAM (>= 8GB recommended for large datasets)
- Required Python packages (see requirements.txt)

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/mdakkasali/Spectral_Flow_Cytometry_Data_Analysis.git
   cd Spectral_Flow_Cytometry_Data_Analysis
   ```

2. Install required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Download the associated dataset from Figshare (link to be provided)

### Reproduction Steps

1. **Environment Setup**: Ensure all dependencies are installed
2. **Data Preparation**: Place downloaded datasets in the appropriate data directory
3. **Sequential Analysis**: Execute Jupyter notebooks in the recommended order:
   - Start with the main analysis notebook
   - Follow the annotated workflow step-by-step
   - Examine intermediate results and quality control outputs
4. **Parameter Customization**: Modify analysis parameters as needed for your specific dataset
5. **Results Generation**: Generate figures and statistical outputs matching the publication

---

## Acknowledgments and Funding

This work was supported by the National Institutes of Health (NIH) through the following grants:

- **National Institute of Allergy and Infectious Diseases (NIAID)**: R01AI134907 (D.J.T.)
- **National Institute of General Medical Sciences (NIGMS)**: T32GM139804 (M.A.A.)
- **National Cancer Institute (NCI)**: P30CA013148 (UAB Comprehensive Cancer Center Core Support)
- **National Institute of Environmental Health Sciences (NIEHS)**: P30ES027797 (UAB Center for Free Radical Biology)

We acknowledge the UAB Flow Cytometry Core Facility, supported by the UAB Comprehensive Cancer Center (P30CA013148), for technical assistance with spectral flow cytometry experiments. We thank the UAB High Performance Computing cluster for computational resources. We also acknowledge the contributions of the UAB Department of Pathology for institutional support and research infrastructure.

Special thanks to all study participants and the research staff who made this work possible. We appreciate the collaborative environment fostered by the UAB Graduate Biomedical Sciences program and the Tyrrell Laboratory members for their valuable input and discussions throughout this project.

---

## Contact

**Md Akkas Ali, PhD Candidate**  
Graduate Biomedical Sciences (GBS)  
Department of Pathology – Molecular & Cellular  
UAB | The University of Alabama at Birmingham  
Office: PBMR2 574 | 901 19th St. South | Birmingham, AL 35205  
P: +1 (205) 381-6106 | mali3@uab.edu

**Connect with me:**
- [Google Scholar](https://scholar.google.com/citations?user=XXXXXXXX)
- [LinkedIn](https://www.linkedin.com/in/md-akkas-ali)
- [X (Twitter)](https://twitter.com/mdakkasali)
- [Personal Website](https://mdakkasali.github.io)
- [Tyrrell Lab](https://www.uab.edu/medicine/tyrrell-lab)

**For questions regarding this repository or publication:**
- **Technical questions**: Create an issue in this GitHub repository
- **Methodological inquiries**: Contact the corresponding author
- **Collaboration opportunities**: Email mali3@uab.edu

We encourage users to report bugs, suggest improvements, and contribute to the ongoing development of this analysis pipeline.

---

## License

This repository supports open science and reproducible research in immunology and flow cytometry analysis. Please refer to the LICENSE file for specific usage terms.

---

*This repository provides comprehensive tools for spectral flow cytometry analysis, enabling researchers to reproduce the methods described in our publication and apply them to their own high-dimensional immune profiling studies.*
