from django.shortcuts import render

# Create your views here.

def index(request):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

    return render(request, 'polarisation/index.html', {})


