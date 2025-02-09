# # Initiate a log
#         timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
#         _fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
#         # Set opath
#         if opath:
#             self.opath = Path(opath)
#             self.opath.mkdir(parents=True, exist_ok=True)
#         else:
#             self.opath = Path.cwd()

#         # Logfile
#         log_file = Path(self.opath / f"lap_Parse_SCNMT_{timestamp}_{uuid.uuid4().hex}.log").name
#         logging.basicConfig(level=logging.INFO)
#         self.logger = logging.getLogger(log_file)
#         file_handler = logging.FileHandler(log_file)
#         # To file
#         file_handler.setFormatter(_fmt)
#         self.logger.addHandler(file_handler)
#         self.logger.propagate = False