"""Minimal DICOM header reader (no pydicom). Reads tags up to PixelData."""
import struct, sys, os, json

LONG_VR = {b'OB', b'OW', b'OF', b'SQ', b'UT', b'UN', b'OD', b'OL', b'UC', b'UR', b'OV'}
NAMES = {
    (0x0008,0x0060):'Modality',(0x0008,0x0070):'Manufacturer',(0x0008,0x1090):'ModelName',
    (0x0008,0x0080):'Institution',(0x0018,0x1000):'DeviceSerial',(0x0018,0x1020):'SoftwareVersions',
    (0x0008,0x0020):'StudyDate',(0x0008,0x0023):'ContentDate',(0x0008,0x0016):'SOPClassUID',
    (0x0008,0x0008):'ImageType',(0x0008,0x0064):'ConversionType',(0x0008,0x103E):'SeriesDescription',
    (0x0008,0x1030):'StudyDescription',(0x0018,0x0015):'BodyPart',(0x0018,0x1030):'ProtocolName',
    (0x0010,0x0010):'PatientName',(0x0010,0x0020):'PatientID',(0x0010,0x0030):'BirthDate',
    (0x0010,0x0040):'Sex',(0x0010,0x1010):'Age',(0x0010,0x1020):'Size',(0x0010,0x1030):'Weight',
    (0x0018,0x0060):'KVP',(0x0018,0x1150):'ExposureTime_ms',(0x0018,0x1151):'TubeCurrent_mA',
    (0x0018,0x1152):'Exposure_mAs',(0x0018,0x1153):'Exposure_uAs',(0x0018,0x8151):'XRayTubeCurrent_uA',
    (0x0018,0x8150):'ExposureTime_us',(0x0018,0x115E):'DAP',(0x0018,0x1405):'RelativeXRayExposure',
    (0x0018,0x1411):'ExposureIndex',(0x0018,0x1413):'TargetExposureIndex',(0x0018,0x1412):'DeviationIndex',
    (0x0018,0x7062):'ExposureControlMode',(0x0018,0x1160):'FilterType',(0x0018,0x7050):'FilterMaterial',
    (0x0018,0x1190):'FocalSpot',(0x0018,0x1110):'SID',(0x0018,0x1111):'SOD',(0x0018,0x1164):'ImagerPixelSpacing',
    (0x0028,0x0030):'PixelSpacing',(0x0018,0x9332):'Exposure_mAs_FD',(0x0018,0x9330):'TubeCurrent_mA_FD',
    (0x0018,0x9328):'ExposureTime_ms_FD',(0x0040,0x0302):'EntranceDose',(0x0040,0x8302):'EntranceDose_mGy',
    (0x0040,0x0316):'OrganDose',(0x0018,0x1508):'PositionerType',(0x0018,0x5101):'ViewPosition',
    (0x0028,0x0010):'Rows',(0x0028,0x0011):'Columns',(0x0028,0x0100):'BitsAllocated',(0x0028,0x0101):'BitsStored',
    (0x0028,0x0004):'Photometric',(0x0002,0x0010):'TransferSyntax',(0x0018,0x1004):'PlateID',
    (0x0018,0x700A):'DetectorID',(0x0018,0x7004):'DetectorType',(0x0018,0x1050):'SpatialResolution',
}
NUMVR = {b'US':'<H',b'UL':'<I',b'SS':'<h',b'SL':'<i',b'FL':'<f',b'FD':'<d'}

def read(path):
    out, priv = {}, {}
    with open(path,'rb') as f:
        data = f.read(400000)
    pos = 132 if data[128:132]==b'DICM' else 0
    explicit = True
    while pos+8 <= len(data):
        g,e = struct.unpack('<HH', data[pos:pos+4])
        if (g,e)==(0x7FE0,0x0010): break
        vr = data[pos+4:pos+6]
        if g!=0x0002 and not explicit:
            vr=None; L=struct.unpack('<I',data[pos+4:pos+8])[0]; hdr=8
        elif vr in LONG_VR:
            L=struct.unpack('<I',data[pos+8:pos+12])[0]; hdr=12
        elif vr.isalpha() and vr.isupper():
            L=struct.unpack('<H',data[pos+6:pos+8])[0]; hdr=8
        else:
            explicit=False; continue
        val_start=pos+hdr
        if L==0xFFFFFFFF:
            # undefined length sequence: skip to sequence delimiter
            d=data.find(b'\xfe\xff\xdd\xe0',val_start)
            pos = d+8 if d>0 else len(data); continue
        raw=data[val_start:val_start+L]
        if vr in NUMVR and L>=struct.calcsize(NUMVR[vr]):
            n=L//struct.calcsize(NUMVR[vr]); val=list(struct.unpack('<'+NUMVR[vr][1]*n,raw))
            val=val[0] if n==1 else val
        elif vr in (b'SQ',b'OB',b'OW',b'UN'):
            val=f'<{vr.decode()} {L}B>'
        else:
            val=raw.decode('latin-1',errors='replace').strip('\x00 ')
        if (g,e)==(0x0002,0x0010):
            explicit = val.strip() != '1.2.840.10008.1.2'
        key=NAMES.get((g,e))
        if key: out[key]=val
        elif g%2==1 and isinstance(val,(str,int,float)) and val!='' :
            priv[f'({g:04X},{e:04X})']=val
        pos=val_start+L
    return out, priv

if __name__=='__main__':
    d=sys.argv[1]; res={}
    for fn in sorted(os.listdir(d), key=lambda x:int(x.split('.')[0])):
        try: res[fn]=read(os.path.join(d,fn))
        except Exception as ex: res[fn]=({'error':str(ex)},{})
    json.dump(res, open(sys.argv[2],'w'), indent=1, default=str)
    print('done', len(res))
