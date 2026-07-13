import numpy as np
import matplotlib.pyplot as plt
from analysis import AnalysisBase
import logging_svs as ls

class AnalysisSatDataStorage(AnalysisBase):
    def __init__(self):
        super().__init__()
        self.ssr_capacity_gbits = 0
        self.initial_fill_gbits = 0
        self.instrument_rate_mbps = 0
        self.downlink_rate_mbps = 0
        self.lat_limit = 90.0
        self.metric = None # Stores [Time, SSR_Level_Gbits, Is_Downlinking, Is_Recording]

    def read_config(self, node):
        self.ssr_capacity_gbits = float(node.find('SSRCapacityGbits').text)
        self.initial_fill_gbits = float(node.find('InitialFillGbits').text)
        self.instrument_rate_mbps = float(node.find('InstrumentRateMbps').text)
        self.downlink_rate_mbps = float(node.find('DownlinkRateMbps').text)
        if node.find('PayloadLatitudeLimit') is not None:
            self.lat_limit = float(node.find('PayloadLatitudeLimit').text)

    def before_loop(self, sm):
        self.metric = np.zeros((sm.num_epoch, 4))
        self.current_fill = self.initial_fill_gbits
        ls.logger.info("Data Storage Analysis Initialized")

    def in_loop(self, sm):
        sat = sm.satellites[0]
        sat.det_lla()
        
        # 1. Recording Logic (matching your Power module latitude trigger)
        is_recording = abs(np.degrees(sat.lla[0])) <= self.lat_limit

        # 2. Downlinking Logic (Ground Station in view)
        is_downlinking = len(sat.idx_stat_in_view) > 0

        # 3. Data Budget Calculation (Mbps * sec / 1000 = Gbits)
        inflow = (self.instrument_rate_mbps / 1000.0) * sm.time_step if is_recording else 0
        outflow = (self.downlink_rate_mbps / 1000.0) * sm.time_step if is_downlinking else 0
        
        self.current_fill += (inflow - outflow)
        
        # Constraints: Cannot be negative, cannot exceed capacity
        self.current_fill = np.clip(self.current_fill, 0, self.ssr_capacity_gbits)

        # Store [DOY, Fill Level, Downlink Status, Record Status]
        self.metric[sm.cnt_epoch, :] = [self.times_f_doy[-1], self.current_fill, 
                                        float(is_downlinking), float(is_recording)]

    def after_loop(self, sm):
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(self.metric[:, 0], self.metric[:, 1], 'g-', label='SSR Fill Level')
        ax1.set_ylabel('Data Stored (Gbits)')
        ax1.set_xlabel('Day of Year (DOY)')
        
        # Shade regions to show activity for visual debugging
        ax1.fill_between(self.metric[:, 0], 0, self.ssr_capacity_gbits, 
                        where=self.metric[:, 2] > 0, color='blue', alpha=0.1, label='Downlink Active')
        
        plt.grid(True)
        plt.legend()
        plt.savefig(sm.output_path('dat_storage.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'ssr_fill_gbits', 'is_downlinking', 'is_recording'],
                       self.metric)

class AnalysisSatDataLatency(AnalysisBase):
    def __init__(self):
        super().__init__()
        self.ssr_capacity_gbits = 0
        self.initial_fill_gbits = 0
        self.instrument_rate_mbps = 0
        self.downlink_rate_mbps = 0
        self.lat_limit = 90.0
        self.ground_proc_min = 0 # New: x minutes for ground processing
        
        self.data_queue = []       
        self.latency_metrics = []  
        self.metric = None         

    def read_config(self, node):
        self.ssr_capacity_gbits = float(node.find('SSRCapacityGbits').text)
        self.initial_fill_gbits = float(node.find('InitialFillGbits').text)
        self.instrument_rate_mbps = float(node.find('InstrumentRateMbps').text)
        self.downlink_rate_mbps = float(node.find('DownlinkRateMbps').text)
        if node.find('PayloadLatitudeLimit') is not None:
            self.lat_limit = float(node.find('PayloadLatitudeLimit').text)
        # Load the ground processing delay (x) from XML
        if node.find('GroundProcessingMin') is not None:
            self.ground_proc_min = float(node.find('GroundProcessingMin').text)

    def before_loop(self, sm):
        self.metric = np.zeros((sm.num_epoch, 4))
        self.current_fill = self.initial_fill_gbits
        if self.initial_fill_gbits > 0:
            self.data_queue.append([sm.time_mjd, self.initial_fill_gbits])
        ls.logger.info(f"Data Analysis Initialized with {self.ground_proc_min}min ground delay")

    def in_loop(self, sm):
        sat = sm.satellites[0]
        sat.det_lla()
        
        is_recording = abs(np.degrees(sat.lla[0])) <= self.lat_limit
        if is_recording:
            generated_gbits = (self.instrument_rate_mbps / 1000.0) * sm.time_step
            self.data_queue.append([sm.time_mjd, generated_gbits])
            self.current_fill += generated_gbits

        is_downlinking = len(sat.idx_stat_in_view) > 0
        if is_downlinking and self.current_fill > 0:
            downlink_capacity = (self.downlink_rate_mbps / 1000.0) * sm.time_step
            
            while downlink_capacity > 0 and len(self.data_queue) > 0:
                packet_time, packet_size = self.data_queue[0]
                
                if packet_size <= downlink_capacity:
                    # Latency = (Time on Orbit) + (Ground Processing x)
                    latency_h = ((sm.time_mjd - packet_time) * 24.0) + (self.ground_proc_min / 60.0)
                    self.latency_metrics.append([self.times_f_doy[-1], latency_h])
                    downlink_capacity -= packet_size
                    self.current_fill -= packet_size
                    self.data_queue.pop(0)
                else:
                    self.data_queue[0][1] -= downlink_capacity
                    self.current_fill -= downlink_capacity
                    latency_h = ((sm.time_mjd - packet_time) * 24.0) + (self.ground_proc_min / 60.0)
                    self.latency_metrics.append([self.times_f_doy[-1], latency_h])
                    downlink_capacity = 0

        self.current_fill = np.clip(self.current_fill, 0, self.ssr_capacity_gbits)
        self.metric[sm.cnt_epoch, :] = [self.times_f_doy[-1], self.current_fill, 
                                        float(is_downlinking), float(is_recording)]

    def after_loop(self, sm):
        self.write_csv(sm, ['doy', 'latency_hours'], self.latency_metrics)
        if not self.latency_metrics:
            return

        lat_data = np.array(self.latency_metrics)
        latencies = lat_data[:, 1] 
        
        # Calculate statistics
        mean_lat = np.mean(latencies)
        p95_lat = np.percentile(latencies, 95)
        max_lat = np.max(latencies)
        
        # Calculate percentage < 2 hours
        pct_under_2h = (np.sum(latencies < 2.0) / len(latencies)) * 100

        fig, (ax2, ax3) = plt.subplots(2, 1, figsize=(12, 16))

        # --- Latency Time Series ---
        ax2.scatter(lat_data[:, 0], latencies, c='red', s=5, alpha=0.3)
        ax2.axhline(2.0, color='black', linestyle=':', label='2h Threshold')
        ax2.set_ylabel('Latency (Hours)')
        ax2.set_title(f'{pct_under_2h:.1f}% of data received in < 2 hours')
        ax2.legend()
        ax2.grid(True)

        # --- Histogram ---
        ax3.hist(latencies, bins=50, color='skyblue', edgecolor='black', alpha=0.7)
        ax3.axvline(mean_lat, color='blue', linestyle='--', label=f'Mean: {mean_lat:.2f}h')
        ax3.axvline(p95_lat, color='orange', linestyle='--', label=f'95%: {p95_lat:.2f}h')
        ax3.set_xlabel('Latency (Hours)')
        ax3.legend()
        ax3.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        plt.savefig(sm.output_path('dat_latency_stats.png'))
        plt.show()