from data.base_dataset import BaseDataset
import os
from pathlib import Path
import uuid
import numpy as np
from utils.read import read_fast5
import torch
from utils.normalization import normalize_signal_from_read_data, med_mad

class BaseFast5Dataset(BaseDataset):
    def __init__(self, opt, manager):
        BaseDataset.__init__(self, opt, manager)
        # Create a Nanopore read from the directory
        opt_data = opt.dataset
        self.recursive = opt_data.recursive
        self.buffer_size = opt_data.buffer_size
        self.window_size = opt_data.window_size
        self.window_overlap = opt_data.window_overlap
        self.trim_signal = opt_data.trim_signal
        if opt_data.fast5_list is None:
            self.data_dir = os.path.join(opt.data_dir, opt.dataset.loc)
            self.data_files = self.find_all_fast5_files()
        else:
            self.data_files = self.read_fast5_list(opt_data.fast5_list)

    def __len__(self):
        return len(self.data_files)

    def __getitem__(self, idx):
        return self.process_reads(self.data_files[idx])

    def find_all_fast5_files(self):
        """Find all fast5 files in a dir recursively
        """
        # find all the files that we have to process
        files_list = list()
        for path in Path(self.data_dir).rglob('*.fast5'):
            files_list.append(str(path))
        files_list = self.buffer_list(files_list, self.buffer_size)
        return files_list

    def read_fast5_list(self, fast5_list):
        """Read a text file with the reads to be processed
        """
        if isinstance(fast5_list, list):
            return self.buffer_list(fast5_list, self.buffer_size)
        # If list store in a text file
        files_list = list()
        with open(fast5_list, 'r') as f:
            for line in f:
                files_list.append(line.strip('\n'))
        files_list = self.buffer_list(files_list, self.buffer_size)
        return files_list

    def buffer_list(self, files_list, buffer_size):
        buffered_list = list()
        for i in range(0, len(files_list), buffer_size):
            buffered_list.append(files_list[i:i+buffer_size])
        # List of list
        return buffered_list

    def normalize(self, read_data):
        return normalize_signal_from_read_data(read_data)

    def trim(self, signal, window_size=40, threshold_factor=2.4, min_elements=3):
        """

        from: https://github.com/nanoporetech/bonito/blob/master/bonito/fast5.py
        """

        min_trim = 10
        signal = signal[min_trim:]

        med, mad = med_mad(signal[-(window_size*100):])

        threshold = med + mad * threshold_factor
        num_windows = len(signal) // window_size

        seen_peak = False

        for pos in range(num_windows):
            start = pos * window_size
            end = start + window_size
            window = signal[start:end]
            if len(window[window > threshold]) > min_elements or seen_peak:
                seen_peak = True
                if window[-1] > threshold:
                    continue
                return min(end + min_trim, len(signal)), len(signal)

        return min_trim, len(signal)

    def chunk(self, signal, chunksize, overlap):
        """
        Convert a read into overlapping chunks before calling

        The first N datapoints will be cut out so that the window ends perfectly
        with the number of datapoints of the read.
        """
        if isinstance(signal, np.ndarray):
            signal = torch.from_numpy(signal)

        T = signal.shape[0]
        if chunksize == 0:
            chunks = signal[None, :]
        elif T < chunksize:
            chunks = torch.nn.functional.pad(signal, (chunksize - T, 0))[None, :]
        else:
            stub = (T - overlap) % (chunksize - overlap)
            chunks = signal[stub:].unfold(0, chunksize, chunksize - overlap)

        return chunks.unsqueeze(1)

    def process_reads(self, read_list):
        """
        Args:
            read_list (list): list of files to be processed

        Returns:
            two arrays, the first one with the normalzized chunked data,
            the second one with the read ids of each chunk.
        """
        chunks_list = list()
        id_list = list()
        l_list = list()

        for read_file in read_list:
            # For every .fast5 file (contains multiple reads)
            reads_data = read_fast5(read_file)
            # A dict: read_id : ReadData Object
            for read_id in reads_data.keys():
                # For every read in a single .fast5 file
                read_data = reads_data[read_id]

                # Return a normalized signal for ReadData.raw : nd.array
                norm_signal = self.normalize(read_data)

                if self.trim_signal:
                    trim, _ = self.trim(norm_signal[:8000])
                    norm_signal = norm_signal[trim:]

                chunks = self.chunk(norm_signal, self.window_size, self.window_overlap)
                # num_chunks * 1 * window_size
                num_chunks = chunks.shape[0]

                # Identify: Which chunk belongs to the same read ?
                uuid_fields = uuid.UUID(read_id).fields
                id_arr = np.zeros((num_chunks, 6), dtype=np.int64)
                for i, uf in enumerate(uuid_fields):
                    id_arr[:, i] = uf

                id_list.append(id_arr)
                l_list.append(np.full((num_chunks,), len(norm_signal)))
                chunks_list.append(chunks)

        out = {
            # x: [n_reads * n_chunks] * window_size
            # id: [n_reads * n_chunks] * 6
            # len: [n_reads * n_chunks]
            'x': torch.vstack(chunks_list).squeeze(1),
            'id': np.vstack(id_list),
            'len': np.concatenate(l_list)
        }
        return out


