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

*[Placeholder for grant/funding agency information to be filled based on publication details]*

We acknowledge the support of [funding agencies and grant numbers], institutional resources, and technical assistance that made this research possible. Detailed acknowledgments can be found in the published manuscript.

---

## Contact

For questions regarding this repository or the associated publication, please:

- **Create an issue** in this GitHub repository for technical questions
- **Contact the corresponding author** for methodological inquiries
- **Email**: [Contact information from publication]

We encourage users to report bugs, suggest improvements, and contribute to the ongoing development of this analysis pipeline.

---

## License

This repository supports open science and reproducible research in immunology and flow cytometry analysis. Please refer to the LICENSE file for specific usage terms.

---

*This repository provides comprehensive tools for spectral flow cytometry analysis, enabling researchers to reproduce the methods described in our publication and apply them to their own high-dimensional immune profiling studies.*
