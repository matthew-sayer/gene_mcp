import os
import dotenv

dotenv.load_dotenv()

NVIDIA_API_KEY = os.environ.get('NVIDIA_API_KEY')
GOOGLE_SEARCH_API_KEY = os.environ.get('GOOGLE_SEARCH_API_KEY')
GOOGLE_CUSTOM_SEARCH_CX = os.environ.get('GOOGLE_CUSTOM_SEARCH_CX')