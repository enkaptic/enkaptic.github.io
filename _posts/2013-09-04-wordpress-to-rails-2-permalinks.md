---
title: WordPress to Rails (2) permalinks
categories:
- Computing
tags:
- Rails
- WordPress
---

Various changes since <a title="WordPress to Rails (1) schema and scaffold" href="/2013/09/03/wordpress-to-rails-1-schema-and-scaffold/">last time</a>, such as adding user login authentication, with the following in the User model, which seems to mimic the WordPress authentication (at least for my password):

{% highlight ruby %}
def check_password(plaintext_password)
  require 'phpass'
  return Phpass.new(8).check(plaintext_password, self.user_pass)
end
{% endhighlight %}

Thanks to <a href="http://blog.amhill.net/2012/03/13/authenticating-against-wordpress-in-rails-3/">Aaron Hill</a> for a helpful post on that topic.

Biggest change has been the introduction of permalinks, and making it possible to view each blog post on a separate page. Within the Post model, I added the following:

{% highlight ruby %}
def self.find_by_permalink(permalink)
  year, month, day, name = permalink.split('/')
  name ||= ''
  find(:first, :conditions => ["YEAR(post_date) = ? AND MONTH(post_date) = ? AND DAYOFMONTH(post_date) = ? AND post_name = ?", year.to_i, month.to_i, day.to_i, name])
end

def to_param
  permalink
end

def permalink
  post_date.strftime("%Y/%m/%d") + (post_name.empty? ? '' : "/#{post_name}")
end
{% endhighlight %}

Thanks to the author of <a href="http://www.dzone.com/snippets/use-contents-wordpress">this snippet</a> for some of that.

Then in <code>routes.rb</code>, I have:

{% highlight ruby %}
resources :posts, :only => [:index, :new, :create]
resources :posts, :only => [:show, :update, :destroy, :edit], :path => '', :id => /\d{4}\/\d{2}\/\d{2}(\/[0-9a-z\-_]+)?/
{% endhighlight %}

This preserves <code>mydomain.com/posts</code> and <code>mydomain.com/posts/new</code> URLs for listing and creating posts, but allows <code>mydomain.com/2013/09/04/my-post</code> as a permalink for the post itself (allowing lower-case letters, numbers, hyphens and underscores).

(I'll leave the smileys in, just to show how clever WordPress is...)

Finally, within <code>PostsController</code>, replace all the <code>Post.find</code> occurrences with <code>Post.find_by_permalink</code>, and I think that's all!