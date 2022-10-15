from .spa import spa as SPA
from .sampa import sampa as SAMPA
from .solpos import solpos as SOLPOS
from .bird import bird as Bird

SPA.__name__ = "SPA"
SAMPA.__name__ = "SAMPA"
SOLPOS.__name__ = "SOLPOS"
Bird.__name__ = "Bird"

__all__ = ["SPA", "SAMPA", "SOLPOS", "Bird"]
