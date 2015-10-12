import sys
import struct
 
if len(sys.argv) < 2:
    usage();
    sys.exit(0)
 
 
def usage():
    print("%s [filename]"%(sys.argv[0]))
 
 
class NotMP4FormatException(Exception):
    pass
 
 
class CHUNK(object):
    def __init__(self, pos):
        self.pos = pos
        self.samples_count = 0
        self.samples_desc_idx = 0
        self.size = 0
 
    def to_annexb(self, filename, sps, pps):
        wf = open(filename, "wb");
        wf.write(struct.pack('>i', 1))
        wf.write(sps[0])
        wf.write(struct.pack('>i', 1))
        wf.write(pps[0])
 
        f = open(sys.argv[1], "rb");
        pos = 0
        f.seek(self.pos)
        while True:
            tsize = struct.unpack('>i', f.read(4))[0]
            if tsize <= 0:
                break
 
            pos += 4
            pos += tsize
 
            if pos > self.size:
                break
 
            data = f.read(tsize)
 
            wf.write(struct.pack('>i', 1))
            wf.write(data)
 
    def __str__(self):
        return "CHUNK(pos: %s, size: %s, samples: %s)"%(self.pos, self.size, self.samples_count)
 
    def __repr__(self):
        return self.__str__()
 
 
class TRACK(object):
    def __init__(self, track_idx, track_atom, mdat_atom):
        self.track_idx = track_idx
        self.track_atom = track_atom
        self.mdat_atom = mdat_atom
        self.stbl_atom = None
        self.ppss = None
        self.spss = None
        self.chunks = []
        self.iframes = None
        self.merge()
 
    def __str__(self):
        return "TRACK(%s)"%(self.track_idx)
 
    def __repr__(self):
        return self.__str__()
 
    def merge(self):
        self.stbl_atom = self.track_atom.find_child_atom("mdia/minf/stbl")
        stco_atom = self.stbl_atom.find_child_atom("stco")
 
        chunks = []
        f = open(sys.argv[1], "rb");
        f.seek(stco_atom.pos+12)
        stco_size = struct.unpack('>i', f.read(4))[0]
 
        for i in range(stco_size):
            p = struct.unpack('>i', f.read(4))[0]
 
            if i != 0:
                before_chunk = chunks[-1]
                before_chunk.size = p - before_chunk.pos
 
            chunks.append(CHUNK(p))
 
        self.chunks = chunks
        stsc_atom = self.stbl_atom.find_child_atom("stsc")
        f.seek(stsc_atom.pos+12)
        stsc_size = struct.unpack('>i', f.read(4))[0]
 
        samples = []
        end_chunk = stco_size
        for i in range(stsc_size):
            start_chunk = struct.unpack('>i', f.read(4))[0]
            sample_count = struct.unpack('>i', f.read(4))[0]
            desc_idx = struct.unpack('>i', f.read(4))[0]
 
            if i != 0:
                samples[i][1] = start_chunk - 1
 
            samples.append([start_chunk, end_chunk+1, sample_count, desc_idx])
 
        idx = 0
        sample_info = samples[idx]
        for i in range(1, stco_size+1):
            chunk = chunks[i-1]
            if i > sample_info[1]:
                idx += 1
                sample_info = samples[idx]
 
            chunk.samples_count = sample_info[2]
            chunk.samples_desc_idx = sample_info[3]
 
        stss_atom = self.stbl_atom.find_child_atom("stss")
        if stss_atom == None:
            return
 
        f.seek(stss_atom.pos+12)
        stss_size = struct.unpack('>i', f.read(4))[0]
 
        iframes = []
        for i in range(stss_size):
            iframes.append(struct.unpack('>i', f.read(4))[0])
 
        self.iframes = iframes
        stsd_atom = self.stbl_atom.find_child_atom("stsd")
        if stsd_atom == None:
            return
 
        self.ppss = stsd_atom.properties[0]["avc"]["pps"]
        self.spss = stsd_atom.properties[0]["avc"]["sps"]
 
 
class ATOM(object):
    def __init__(self, size, name, pos):
        self.size = size
        self.name = name
        self.pos = pos
        self.children = []
        self.properties = None
 
    def find_child_atom_internal(self, atoms, part_arr):
        name = part_arr[0]
        for atom in atoms:
            if atom.name == name:
                if len(part_arr) == 1:
                    return atom
 
                return self.find_child_atom_internal(atom.children, part_arr[1:])
 
        return None
 
    def find_child_atom(self, name):
        part_arr = name.split("/")
        return self.find_child_atom_internal(self.children, part_arr)
 
    def __str__(self):
        return "%s(%s)"%(self.name, self.size)
 
    def __repr__(self):
        return self.__str__()
 
 
class MP4(object):
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename, "rb");
        self.children = []
        self.track_size = 0
        self.moov_atom = None
        self.mdat_atom = None
        self.tracks = []
 
    def is_parent_atom(self, name):
        return name not in ['mdat', 'tkhd', 'vmhd']
 
    def create_empty_atom(self):
        return ATOM(0, "", 0)
 
    def get_moov_atom(self):
        for atom in self.children:
            if atom.name == "moov":
                return atom
 
        raise NotMP4FormatException()
 
    def get_mdat_atom(self):
        for atom in self.children:
            if atom.name == "mdat":
                return atom
 
        raise NotMP4FormatException()
 
    def get_track_size(self):
        return self.track_size
 
    def get_track_size_internal(self, atom):
        count = 0
        for atom in atom.children:
            if atom.name == "trak":
                count += 1
 
        return count;
 
    def parse(self, start_pos = 0):
        #mp4 container follow BIG ENDIAN
        next_pos = start_pos
 
        try:
            while True:
                atom = self.parse_internal(next_pos)
                self.children.append(atom)
                next_pos += atom.size
 
        except struct.error:
            pass
 
        except:
            raise NotMP4FormatException()
 
        self.moov_atom = self.get_moov_atom()
        self.mdat_atom = self.get_mdat_atom()
        self.track_size = self.get_track_size_internal(self.moov_atom)
 
        self.tracks = self.merge_tracks()
        return True
 
    def merge_tracks(self):
        tracks = []
        count = 0
        for atom in self.moov_atom.children:
            if atom.name == "trak":
                tracks.append(TRACK(count, atom, self.mdat_atom))
                count += 1
 
        return tracks
 
    def traverse(self, udf = None):
        self.traverse_internal(self.children, 0, udf)
 
    def traverse_internal(self, atoms, depth, udf = None):
        buf = ""
        for i in range(depth):
            buf += "    "
 
        for atom in atoms:
            print "%s%s"%(buf, atom)
            if udf is not None:
                udf(atom)
 
            self.traverse_internal(atom.children, depth+1, udf)
 
    def get_atom(self, pos):
        self.f.seek(pos)
        size = struct.unpack('>i', self.f.read(4))[0]
        name = self.f.read(4)
        return ATOM(size, name, pos)
 
    def parse_avcC(self, avc, name, size):
        avcC = {}
        spss = []
        ppss = []
        version = struct.unpack('>b', self.f.read(1))[0]
        avc_profile_idc = struct.unpack('>b', self.f.read(1))[0]
        profile_compatibility = struct.unpack('>b', self.f.read(1))[0]
        avc_level_idc = struct.unpack('>b', self.f.read(1))[0]
 
        lengh_size_minus_one = (struct.unpack('>b', self.f.read(1))[0]) & 0x03 + 1
        num_of_sps = (struct.unpack('>b', self.f.read(1))[0]) & 0x1F
        for i in range(num_of_sps):
            length_sps = struct.unpack('>h', self.f.read(2))[0]
            sps = self.f.read(length_sps)
            spss.append(sps)
 
        num_of_pps = struct.unpack('>b', self.f.read(1))[0]
        for i in range(num_of_pps):
            length_pps = struct.unpack('>h', self.f.read(2))[0]
            pps = self.f.read(length_pps)
            ppss.append(pps)
        
        avcC["length_size_minus_one"] = lengh_size_minus_one
        avcC["sps"] = spss
        avcC["pps"] = ppss
        return avcC
 
    def parse_avc_internal(self, atom):
        avc = {}
        size = struct.unpack('>i', self.f.read(4))[0]
        name = self.f.read(4)
        if name != "avc1":
            return None
 
        avc["name"] = name
        self.f.read(24)
        avc["w"] = struct.unpack('>h', self.f.read(2))[0]
        avc["h"] = struct.unpack('>h', self.f.read(2))[0]
        avc["hres"] = struct.unpack('>i', self.f.read(4))[0]
        avc["vres"] = struct.unpack('>i', self.f.read(4))[0]
        self.f.read(4)
 
        frame_count = struct.unpack('>h', self.f.read(2))[0]
        if frame_count != 1:
            return None
 
        self.f.read(32)
        depth = struct.unpack('>h', self.f.read(2))[0]
        if depth != 0x18:
            return None
 
        pd = struct.unpack('>h', self.f.read(2))[0]
        if pd != -1:
            return None
 
        while True:
            tsize = struct.unpack('>i', self.f.read(4))[0]
            tname = self.f.read(4)
 
            if tname == "avcC":
                avc["avc"] = self.parse_avcC(avc, tname, tsize)
                break
            else:
                self.f.read(tsize-8)
        
        return avc
 
    def parse_avc(self, atom):
        self.f.seek(atom.pos+12)
        entry_count = struct.unpack('>i', self.f.read(4))[0]
        entries = []
 
        for i in range(entry_count):
            entry = self.parse_avc_internal(atom)
            if entry is not None:
                entries.append(entry)
 
        return entries
        
 
    def parse_internal(self, pos, total_size = 0):
        atom = self.get_atom(pos)
        if total_size > 0 and atom.size > total_size:
            return self.create_empty_atom()
 
        if self.is_parent_atom(atom.name) == False:
            return atom
 
        if atom.name == "stsd":
            child = self.parse_avc(atom)
            atom.properties = child
            return atom
 
        next_pos = atom.pos + 8
        temp_size = atom.size
 
        while (next_pos+8) < (atom.pos + atom.size):
            child = self.parse_internal(next_pos, atom.size)
            if (child.size >= atom.size) or child.size <= 0:
                break
 
            atom.children.append(child)
            next_pos += child.size
 
        return atom
 
 
def buffer_to_lines(buf):
    arr = []
    size = len(buf)
    mod = size % 16
    line = ""
    for i in range(size):
        if i != 0 and i % 16 == 0:
            arr.append(line)
            line = ""
        else:
            line += " "
 
        val = ord(buf[i])
        v = "{0:02x}".format(val)
        line += v
 
    if mod != 0:
        for i in range(16-mod):
            line += " 00"
 
    if len(line) != 0:
        arr.append(line)
 
    return arr
 
 
if __name__ == "__main__":
    mp4 = MP4(sys.argv[1])
    mp4.parse()
    mp4.traverse()
 
    track = mp4.tracks[0]
#   print track.chunks
    track.chunks[0].to_annexb("/Users/charsyam/chunk0-0", track.spss, track.ppss)
#   track.chunks[1].to_annexb("d:\\chunk1", track.spss, track.ppss, size)
