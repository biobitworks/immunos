---
source: /Users/byron/projects/data/immunos_data/research/scifact/script/download-data.sh
relative: data/immunos_data/research/scifact/script/download-data.sh
generated_at: 2025-12-23 10:28
---

```bash
#!/bin/bash

# Check if the data directory already exists
if [ -e "data" ]
then
    echo "Data directory already exists. Skip download."
    exit 0
fi

# Download and unpack data from AI2 S3 bucket.
name="https://scifact.s3-us-west-2.amazonaws.com/release/latest/data.tar.gz"
wget $name

tar -xvf data.tar.gz
rm data.tar.gz

```
