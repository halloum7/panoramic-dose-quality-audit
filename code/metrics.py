"""Objective no-reference image-quality metrics for panoramic DICOM (8-bit, uncompressed)."""
import os, sys, struct, json
import numpy as np
from scipy import ndimage
from skimage.measure import shannon_entropy

def load_pixels(path):
    with open(path, 'rb') as f:
        data = f.read()
    # find PixelData tag (7FE0,0010) explicit VR OB/OW
    i = data.rfind(b'\xe0\x7f\x10\x00')
    vr = data[i+4:i+6]
    if vr in (b'OB', b'OW'):
        L = struct.unpack('<I', data[i+8:i+12])[0]; start = i + 12
    else:
        L = struct.unpack('<I', data[i+4:i+8])[0]; start = i + 8
    # rows / cols
    def us(tag):
        j = data.find(tag)
        return struct.unpack('<H', data[j+8:j+10])[0]
    rows = us(b'\x28\x00\x10\x00US'); cols = us(b'\x28\x00\x11\x00US')
    px = np.frombuffer(data[start:start+rows*cols], dtype=np.uint8).reshape(rows, cols)
    return px

def immerkaer_sigma(img):
    """Fast noise std estimate (Immerkaer 1996)."""
    M = np.array([[1, -2, 1], [-2, 4, -2], [1, -2, 1]], dtype=float)
    H, W = img.shape
    conv = ndimage.convolve(img.astype(float), M)
    s = np.sum(np.abs(conv[1:-1, 1:-1]))
    return s * np.sqrt(0.5 * np.pi) / (6.0 * (W - 2) * (H - 2))

def metrics(px):
    H, W = px.shape
    roi = px[int(0.15*H):int(0.85*H), int(0.10*W):int(0.90*W)].astype(float)  # dento-alveolar region
    mean = roi.mean()
    p5, p95 = np.percentile(roi, [5, 95])
    sigma = immerkaer_sigma(roi)
    contrast = p95 - p5
    smooth = ndimage.gaussian_filter(roi, 1.0)
    gx = ndimage.sobel(smooth, 1); gy = ndimage.sobel(smooth, 0)
    tenengrad = np.mean(gx**2 + gy**2)
    return dict(
        mean=round(mean, 1),
        noise_sigma=round(sigma, 2),
        contrast_p95_p5=round(contrast, 1),
        snr=round(mean / sigma, 1),
        cnr=round(contrast / sigma, 1),
        sharpness_tenengrad=round(tenengrad / 1000, 2),
        entropy=round(shannon_entropy(roi.astype(np.uint8)), 3),
        sat_white_pct=round(100 * np.mean(px >= 250), 2),
        sat_black_pct=round(100 * np.mean(px <= 5), 2),
    )

if __name__ == '__main__':
    d, files_json, out = sys.argv[1], sys.argv[2], sys.argv[3]
    files = json.load(open(files_json))
    res = {}
    for fn in files:
        res[fn] = metrics(load_pixels(os.path.join(d, fn)))
    json.dump(res, open(out, 'w'), indent=1)
    print('done', len(res))
