""" 
    This script, adapted from trufor_test.py, integrates TruFor into Sherloq for straightforward use. 
    TruFor is an AI-driven solution designed for Digital Image Forensics. Although powerful, AI approaches 
    may not always be as reliable as their performance statistics might suggest. Nonetheless, TruFor can assist 
    forensic analysts and provide evidence regarding the authenticity of digital images.

    Original TruFor Work:
    Research Group of University Federico II of Naples ('GRIP-UNINA')
    https://github.com/grip-unina/TruFor

    Reference Bibtex:
    @InProceedings{Guillaro_2023_CVPR,
        author    = {Guillaro, Fabrizio and Cozzolino, Davide and Sud, Avneesh and Dufour, Nicholas and Verdoliva, Luisa},
        title     = {TruFor: Leveraging All-Round Clues for Trustworthy Image Forgery Detection and Localization},
        booktitle = {Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)},
        month     = {June},
        year      = {2023},
        pages     = {20606-20615}
    }

    Created in September 2024
    @author: github user
"""

import sys, os
import argparse
import numpy as np
from tqdm import tqdm
from glob import glob

import torch
from torch.nn import functional as F

path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
if path not in sys.path:
    sys.path.insert(0, path)

from TruFor_main.test_docker.src.config import update_config
from TruFor_main.test_docker.src.config import _C as config
from TruFor_main.test_docker.src.data_core import myDataset

def process_image(input_path, gpu):

    parser = argparse.ArgumentParser(description='Test TruFor')
    parser.add_argument('-gpu', '--gpu', type=int, default=0, help='device, use -1 for cpu')
    parser.add_argument('-in', '--input', type=str, default='../images',
                        help='can be a single file, a directory or a glob statement')
    parser.add_argument('-out', '--output', type=str, default='../output', help='output folder')
    parser.add_argument('-save_np', '--save_np', action='store_true', help='whether to save the Noiseprint++ or not')
    parser.add_argument('opts', help="other options", default=None, nargs=argparse.REMAINDER)

    args = parser.parse_args()
    update_config(config, args)


    # Set device (GPU or CPU)
    device = 'cuda:%d' % gpu if gpu >= 0 else 'cpu'
    np.set_printoptions(formatter={'float': '{: 7.3f}'.format})

    if device != 'cpu':
        # cudnn setting
        import torch.backends.cudnn as cudnn

        cudnn.benchmark = config.CUDNN.BENCHMARK
        cudnn.deterministic = config.CUDNN.DETERMINISTIC
        cudnn.enabled = config.CUDNN.ENABLED

    # Resolve input (file or folder)
    if '*' in input_path:
        list_img = glob(input_path, recursive=True)
        list_img = [img for img in list_img if not os.path.isdir(img)]
    elif os.path.isfile(input_path):
        list_img = [input_path]
    elif os.path.isdir(input_path):
        list_img = glob(os.path.join(input_path, '**/*'), recursive=True)
        list_img = [img for img in list_img if not os.path.isdir(img)]
    else:
        raise ValueError("Input is neither a file nor a folder")

    test_dataset = myDataset(list_img=list_img)

    testloader = torch.utils.data.DataLoader(
        test_dataset,
        batch_size=1  # Batch size of 1 to allow arbitrary input sizes
    )

    # Load model
    if config.TEST.MODEL_FILE:
        model_state_file = config.TEST.MODEL_FILE
    else:
        raise ValueError("Model file is not specified.")

    print(f'=> loading model from {model_state_file}')
    checkpoint = torch.load(model_state_file, map_location=torch.device(device))

    if config.MODEL.NAME == 'detconfcmx':
        from models.cmx.builder_np_conf import myEncoderDecoder as confcmx
        model = confcmx(cfg=config)
    else:
        raise NotImplementedError('Model not implemented')

    model.load_state_dict(checkpoint['state_dict'])
    model = model.to(device)

    # Calculate output
    with torch.no_grad():
        for index, (rgb, path) in enumerate(tqdm(testloader)):
            path = path[0]

            try:
                rgb = rgb.to(device)
                model.eval()

                pred, conf, det, npp = model(rgb)

                det_score = torch.sigmoid(det).item()

                pred = torch.squeeze(pred, 0)
                pred = F.softmax(pred, dim=0)[1]
                pred = pred.cpu().numpy()
                torch.cuda.empty_cache()#manage CUDA memory usage by freeing GPU memory after use 
                # Return the prediction map
                return pred, det_score

            except Exception as e:
                print(f"Exception during prediction: {e}")
                pass

    return None, None  # Return None if no prediction was made
