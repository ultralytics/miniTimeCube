<a href="https://www.ultralytics.com/"><img src="https://raw.githubusercontent.com/ultralytics/assets/main/logo/Ultralytics_Logotype_Original.svg" width="320" alt="Ultralytics logo"></a>

# 🌟 Introduction

Welcome to the official repository for the miniTimeCube (mTC) simulation and analysis code. This documentation guides you through the setup process, describes how to run simulations, and provides an overview of the analysis capabilities included in this project. Whether you are a physicist, a [data scientist](https://www.ultralytics.com/glossary/data-analytics), or an enthusiast in the field of particle detection, you'll find the tools and information needed to work with the mTC.

[![Ultralytics Actions](https://github.com/ultralytics/miniTimeCube/actions/workflows/format.yml/badge.svg)](https://github.com/ultralytics/miniTimeCube/actions/workflows/format.yml)
[![Ultralytics Discord](https://img.shields.io/discord/1089800235347353640?logo=discord&logoColor=white&label=Discord&color=blue)](https://discord.com/invite/ultralytics)
[![Ultralytics Forums](https://img.shields.io/discourse/users?server=https%3A%2F%2Fcommunity.ultralytics.com&logo=discourse&label=Forums&color=blue)](https://community.ultralytics.com/)
[![Ultralytics Reddit](https://img.shields.io/reddit/subreddit-subscribers/ultralytics?style=flat&logo=reddit&logoColor=white&label=Reddit&color=blue)](https://reddit.com/r/ultralytics)

## 📚 Description

The miniTimeCube (mTC) is an innovative, compact neutrino detector designed as part of a larger scientific endeavor to study and understand neutrino interactions. This repository contains all the necessary simulation and analysis code related to the mTC project.

The scientific paper associated with this project, titled **"Invited Article: miniTimeCube,"** was authored by a team of respected researchers. It presents a deep dive into the design, operation, and capabilities of the mTC. You can access the full article for a comprehensive understanding via its [DOI link](http://dx.doi.org/10.1063/1.4942243). Understanding concepts like [object detection](https://docs.ultralytics.com/tasks/detect/) can be helpful when working with particle detection data.

<div align="center">
  <img src="https://github.com/ultralytics/mtc/blob/main/cover.jpg" alt="mTC">
</div>

## 🛠 Requirements

To work with the miniTimeCube simulation and analysis code, you will need:

-   **MATLAB** version 2018a or newer. MATLAB is a high-level language and interactive environment for numerical computation, visualization, and programming. If you do not have MATLAB installed, please visit the [official MATLAB website](https://www.mathworks.com/products/matlab.html). Familiarity with [data visualization](https://www.ultralytics.com/glossary/data-visualization) techniques is also beneficial.

In addition, please make sure to clone the common functions repository and add it to your MATLAB path:

```bash
# Clone the functions-matlab repository using Git
git clone https://github.com/ultralytics/functions-matlab
```

```matlab
% Add the repository to the MATLAB path
addpath(genpath('/path/to/functions-matlab')) % Replace /path/to/ with the actual path
```

Ensure that the following toolboxes are installed in your MATLAB environment:

-   **Statistics and Machine Learning Toolbox**: Provides functions and apps to describe, analyze, and model data using statistics and [machine learning](https://www.ultralytics.com/glossary/machine-learning-ml).
-   **Signal Processing Toolbox**: Offers a variety of tools and algorithms for signal processing tasks, crucial for analyzing detector outputs.

## 🏃 Running the Code

To execute the mTC simulation and analysis code within MATLAB, simply enter the following command in your MATLAB terminal:

```matlab
% Launch the nView interface
nView
```

This command starts the `nView` interface, where you can interact with the simulation and analysis tools provided for the miniTimeCube. For insights into optimizing performance, consider exploring resources on [real-time inference](https://www.ultralytics.com/glossary/real-time-inference).

## 🤝 Contribute

We welcome contributions from the community! Whether you're fixing bugs, adding new features, or improving documentation, your input is invaluable. Take a look at our [Contributing Guide](https://docs.ultralytics.com/help/contributing/) to get started. Also, we'd love to hear about your experience with Ultralytics products. Please consider filling out our [Survey](https://www.ultralytics.com/survey?utm_source=github&utm_medium=social&utm_campaign=Survey). A huge 🙏 and thank you to all of our contributors!

[![Ultralytics open-source contributors](https://raw.githubusercontent.com/ultralytics/assets/main/im/image-contributors.png)](https://github.com/ultralytics/ultralytics/graphs/contributors)

## ©️ License

Ultralytics offers two licensing options to accommodate diverse needs:

-   **AGPL-3.0 License**: Ideal for students and enthusiasts, this [OSI-approved](https://opensource.org/license/agpl-3-0/) open-source license promotes collaboration and knowledge sharing. See the [LICENSE](https://github.com/ultralytics/miniTimeCube/blob/main/LICENSE) file for details.
-   **Enterprise License**: Designed for commercial use, this license permits seamless integration of Ultralytics software and AI models into commercial products and services. For commercial applications, please contact us via [Ultralytics Licensing](https://www.ultralytics.com/license).

## 📬 Contact Us

For bug reports, feature requests, and contributions, please visit [GitHub Issues](https://github.com/ultralytics/miniTimeCube/issues). For questions, discussions, and community support regarding this project and other Ultralytics initiatives, join our [Discord](https://discord.com/invite/ultralytics)!

<br>
<div align="center">
  <a href="https://github.com/ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-github.png" width="3%" alt="Ultralytics GitHub"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://www.linkedin.com/company/ultralytics/"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-linkedin.png" width="3%" alt="Ultralytics LinkedIn"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://twitter.com/ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-twitter.png" width="3%" alt="Ultralytics Twitter"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://youtube.com/ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-youtube.png" width="3%" alt="Ultralytics YouTube"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://www.tiktok.com/@ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-tiktok.png" width="3%" alt="Ultralytics TikTok"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://ultralytics.com/bilibili"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-bilibili.png" width="3%" alt="Ultralytics BiliBili"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://discord.com/invite/ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-discord.png" width="3%" alt="Ultralytics Discord"></a>
</div>
