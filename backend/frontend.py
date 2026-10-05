"""Serve known React entry routes without masking API or asset errors."""
from fastapi.staticfiles import StaticFiles


class FrontendFiles(StaticFiles):
    async def get_response(self, path, scope):
        if path.strip('/') in {'about', 'species', 'resources', 'tools'}:
            path = 'index.html'
        return await super().get_response(path, scope)
