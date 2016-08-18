---
title: bundle open gem source
date: 2015-01-12 12:41:24.966441000 +00:00
---
I didn't realise you could do this, from a Ruby or Rails project using [bundler](http://bundler.io/)/<code>Gemfile</code>:

{% highlight bash %}
bundle open activerecord
{% endhighlight %}

Then it opens the source code for that gem (<code>activerecord</code> in this case), using the version specified in your <code>Gemfile.lock</code>.

Thanks to [Justin Weiss](http://www.justinweiss.com/) for this tip, from his weekly email.

You can specify your favourite editor in your <code>~/.bash_profile</code> (or similar). E.g., for [Sublime Text](https://www.sublimetext.com/):

{% highlight bash %}
export EDITOR='subl -w'
{% endhighlight %}

(NB: <code>-w</code> makes it wait for the files to be closed before returning. This is worth doing, as some commands will pause until control returns, and then continue, for example, <code>git commit</code>.)
