# Use the Python 3.12 slim image to save some space
FROM python:3.12-slim

# Stop Python from generating .pyc files
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install the Linux system, GUI, and Chromium dependencies to satisfy the requirements of PySide6
RUN apt-get update && apt-get install -y \
    libgl1 \
    libglib2.0-0 \
    libxkbcommon-x11-0 \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-keysyms1 \
    libxcb-randr0 \
    libxcb-render-util0 \
    libxcb-shape0 \
    libxcb-xinerama0 \
    libxcb-xkb1 \
    libx11-xcb1 \
    libfontconfig1 \
    libdbus-1-3 \
    libmagic-dev \
    libegl1 \
    libxcb-cursor0 \
    libnss3 \
    libxcomposite1 \
    libxdamage1 \
    libxrandr2 \
    libasound2 \
    libxtst6 \
    libxi6 \
    libxfixes3 \
    libxss1 \
    libxkbfile1 \
    libgssapi-krb5-2 \
    libatk1.0-0 \
    libatk-bridge2.0-0 \
    libcups2 \
    libdrm2 \
    libgbm1 \
    libpango-1.0-0 \
    libcairo2 \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory to the root app folder first to copy the requirements files
WORKDIR /app

# Copy both requirements files from the gui folder
COPY gui/requirements.txt gui/requirements_ai_solutions.txt ./

# Install setuptools, standard requirements, AI requirements, and force-upgrade XGBoost to ensure compatibility
RUN pip install --no-cache-dir setuptools && \
    pip install --no-cache-dir -r requirements.txt -r requirements_ai_solutions.txt \
    --extra-index-url https://download.pytorch.org/whl/cpu && \
    pip install --no-cache-dir --upgrade xgboost

# Copy everything else (the Sherloq app) into /app
COPY . .

# Change the working directory to gui so the relative icon paths work correctly
WORKDIR /app/gui

# Start the application
CMD ["python", "sherloq.py"]