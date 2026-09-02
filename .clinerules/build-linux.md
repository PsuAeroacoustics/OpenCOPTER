---
description: Workflow for building OpenCOPTER for debug on linux
author: Rob Rau
version: 1.0
category: "Cline Core"
tags: ["build", "process", "core-behavior"]
globs: ["*"]
---

I MUST try to build the project to ensure there are no compilation errors. I MUST use the following steps to ensure the build succeeds:
* I MUST locate `conda`
* I MUST run `conda activate opencopter` in the root of the repository
* I MUST run `./build_linux.sh native debug` in the root of the repository

These steps ensure that the entirety of OpenCOPTER is built successfully in debug mode.
